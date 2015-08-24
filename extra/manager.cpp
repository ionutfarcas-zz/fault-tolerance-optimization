#include "combigrid/utils/LevelVector.hpp"
#include "combigrid/utils/Types.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"
#include "boost/lexical_cast.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include "combigrid/loadmodel/AnisotropyLoadModel.hpp"
#include "combigrid/loadmodel/LinearLoadModel.hpp"
#include "combigrid/mpi/MPISystem.hpp"
#include "combigrid/manager/TaskLog.hpp"
#include "combigrid/manager/ProcessGroupManager.hpp"
#include "combigrid/task/GeneTask.hpp"
#include "combigrid/sparsegrid/SGrid.hpp"
#include "combigrid/fullgrid/FullGrid.hpp"
#include "combigrid/manager/ProcessManager.hpp"
#include "combigrid/manager/ProcessGroupSignals.hpp"

#include "opt_combi_technique/inc/opt_combi_technique.hpp"

using namespace combigrid;

#define SPCFILENAME "spaces.dat"
#define NMAXFILENAME "nmax.dat"

#ifdef DEBUG
volatile int endStallForDebugger=0;
#else
volatile int endStallForDebugger = 1;
#endif

void
writeGroupTimes( const char* filename, ProcessGroupManagerContainer& pgroups );
void
writeOverallTime( const char* filename, real time, size_t ngroups,
  size_t ntasks );
void
writeSolutionTimes( const char* filename, GeneTaskContainer& tasks );
void
writeTaskTimes( const char* outputdir, GeneTaskContainer& tasks );

void writeEVEvolution( const char* filename_gamma, const char* filename_omega,
 GeneTaskContainer& tasks );

void
stallForDebugger()
{
  while (!endStallForDebugger)
    ;
}

bool sortTaskByID( const GeneTask &t1, const GeneTask &t2 ){
  return ( t1.getID() < t2.getID() );
}

int
main( int argc, char** argv )
{
  MPI_Init( &argc, &argv );

  stallForDebugger();

  // read parameter file
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini( "ctparam", cfg );

  size_t ngroup = cfg.get<size_t>( "manager.ngroup" );
  size_t nprocs = cfg.get<size_t>( "manager.nprocs" );
  std::string basename = cfg.get<std::string>( "manager.basename" );
  size_t ntimesteps_ev_calc = cfg.get<size_t>( "ct.ntimesteps_ev_calc" );

  std::string fg_file_prefix = cfg.get<std::string>( "ct.fg_file_prefix" );
  real dt = cfg.get<real>( "ct.dt_max" );

  NormalizationType normalization_type = cfg.get<NormalizationType>(
    "ct.normalization_type" );
  int normalization_p = cfg.get<int>( "ct.normalization_p" );

  // create local communicator
  int grank, gsize;
  MPI_Comm_rank( MPI_COMM_WORLD, &grank );
  MPI_Comm_size( MPI_COMM_WORLD, &gsize );
  int managerID = gsize - 1;
  MPI_Comm lcomm, pcomm;
  int color = grank / nprocs;
  int key = grank - color * nprocs;
  MPI_Comm_split( MPI_COMM_WORLD, color, key, &lcomm );

  // Gene creates another comm which we do not need, but it is necessary
  // to execute comm_split again
  MPI_Comm_split( MPI_COMM_WORLD, key, color, &pcomm );

  // create global communicator containing manager and pgroup roots
  MPI_Group worldGroup;
  MPI_Comm_group( MPI_COMM_WORLD, &worldGroup );

  std::vector<int> ranks( ngroup + 1 );
  for( size_t i = 0; i < ngroup; ++i )
    ranks[i] = i * nprocs;
  ranks.back() = managerID;

  MPI_Group rootGroup;
  MPI_Group_incl( worldGroup, int( ranks.size() ), &ranks[0], &rootGroup );

  CommunicatorType gcomm;
  MPI_Comm_create( MPI_COMM_WORLD, rootGroup, &gcomm );
  int newsize;
  MPI_Comm_size( gcomm, &newsize );
  managerID = newsize - 1;

  int rold, rnew;
  MPI_Comm_rank( MPI_COMM_WORLD, &rold );
  MPI_Comm_rank( gcomm, &rnew );
  std::cout << "manager: rank old " << rold << " rnew " << rnew << " managerID "
  << managerID << std::endl;

  // create ProcessGroups
  ProcessGroupManagerContainer pgroups;
  for( size_t i = 0; i < ngroup; ++i ){
    // todo: order of ranks in new group?
    int pgroupRootID = static_cast<int>( i );
    ProcessGroupManager grp( managerID, pgroupRootID, gcomm );
    pgroups.push_back( grp );
  }

  // create load model
  LoadModel* loadmodel = new AnisotropyLoadModel();

  // read nmax and lmin from file
  std::ifstream nmaxfile( NMAXFILENAME );
  std::string line;
  LevelVector nmax( 5 );
  LevelVector lmin( 5 );
  LevelVector leval( 5 );

  size_t ntimesteps_total;
  size_t ntimesteps_combi;
  {
    std::getline( nmaxfile, line );
    std::stringstream ss( line );
    ss >> nmax[0] >> nmax[1] >> nmax[2] >> nmax[3] >> nmax[4];
  }
  {
    std::getline( nmaxfile, line );
    std::stringstream ss( line );
    ss >> lmin[0] >> lmin[1] >> lmin[2] >> lmin[3] >> lmin[4];
  }
  {
    std::getline( nmaxfile, line );
    std::stringstream ss( line );
    ss >> leval[0] >> leval[1] >> leval[2] >> leval[3] >> leval[4];
  }
  {
    std::getline( nmaxfile, line );
    std::stringstream ss( line );
    ss >> ntimesteps_total >> ntimesteps_combi;
  }

  assert( ntimesteps_total >= ntimesteps_combi );
  assert( ntimesteps_total % ntimesteps_combi == 0 );
  if( ntimesteps_combi > ntimesteps_ev_calc )
    assert( ntimesteps_combi % ntimesteps_ev_calc == 0 );
  else
    assert( ntimesteps_ev_calc % ntimesteps_combi == 0 );

  // set boundary vector
  std::vector<bool> boundary( 5 );
  boundary[0] = true;
  boundary[1] = false;
  boundary[2] = true;
  boundary[3] = false;
  boundary[4] = true;

  // determine number of time steps of one run
  // todo: find a better solution here
  size_t nsteps =
  ( ntimesteps_combi <= ntimesteps_ev_calc ) ?
  ntimesteps_combi : ntimesteps_ev_calc;

  // create GeneInstances - read levelvectors from file
  GeneTaskContainer tasks;
  std::ifstream spcfile( SPCFILENAME );

  /* optimization code related parameters */
  double level_coeff = 0.0;

  vec2d levels;
  std::vector<int> l_max(lmin.begin(), lmin.end());
  std::vector<int> l_min(nmax.begin(), nmax.end());

  combi_grid_dict given_dict;

  levels.push_back(l_min);
  levels.push_back(l_max);

  vec2d faults {{2, 1, 1, 1, 1}};
  /****************************************/

  while (std::getline( spcfile, line ))
  {
    std::stringstream ss( line );
    int id;
    LevelVector l( 5 );
    combigrid::real coeff;
    ss >> id >> l[0] >> l[1] >> l[2] >> l[3] >> l[4] >> coeff;

    /********** opt code related ***********/

    std::vector<int> levels(l.begin(), l.end());
    level_coeff = static_cast<double>(coeff);
    given_dict.insert(std::make_pair(levels, level_coeff));

    /***************************************/

    // path to gene instance
    std::stringstream ss2;
    ss2 << "../" << basename << id;
    std::string path = ss2.str();

    GeneTask t( l, nmax, lmin, boundary, coeff, path, loadmodel, dt, nsteps );
    tasks.push_back( t );
  }
  spcfile.close();

  // create Manager with process groups
  ProcessManager manager( pgroups, tasks, managerID, gcomm );
  
  double tstart = MPI_Wtime();  

  std::vector<int> leval_tmp2( leval.begin(), leval.end() );
  for( size_t i = 0; i < ntimesteps_total; )
  {
    if( i == 0 )
    {
      // init simulations and run for first time steps
      manager.runfirst();
    }
    else
    {
      manager.runnext();
    }

    i += nsteps;

    if( i % ntimesteps_combi == 0 )
    {
      /*************************************/  
      lp_opt::LP_OPT_INTERP opt_interp(
        levels, 
        5,
        GLP_MAX,
        given_dict,
        faults);
      /*************************************/

      FullGrid<complex> fg( 5, leval_tmp2, boundary );
      manager.combineFG( fg, normalization_type, normalization_p );
    }

    if( i % ntimesteps_ev_calc == 0 )
    {
      FullGrid<complex> fg( 5, leval_tmp2, boundary );
      manager.gridEval( fg );

      // write fg solution
      std::string filename = fg_file_prefix
      + boost::lexical_cast<std::string>( i ) + ".dat";
      fg.save( filename );

      std::cout << "lp norm of fg after grideval"
      << fg.getlpNorm( normalization_p ) << std::endl;
    }
  }

  double time = MPI_Wtime() - tstart;

  // sync tasks
  manager.sync();

  // send exit signal to all process groups
  manager.exit();

  // the order in the tasks container may change due to scheduling
  // for output we sort it by task ID, so that it's consistent with spaces.dat
  std::sort( tasks.begin(), tasks.end(), sortTaskByID );

  /*
  manager.writeCombinationTimes( "../out/combinationTimes.dat" );
  manager.writeSyncTimes( "../out/syncTimes.dat" );
  manager.writeGridEvalTimes( "../out/gridEvalTimes.dat" );
  */

  std::cout << "overall runtime = " << time << std::endl;

  // write times of process groups to file
  //writeGroupTimes( "../out/groupTimes.dat", pgroups );

  // write overall time to file
  //writeOverallTime( "../out/overallTime.dat", time, pgroups.size(), tasks.size() );

  // write solution times to file
  //writeSolutionTimes( "../out/solutionTimes.dat", tasks );

  // write runtimes of each task
  //writeTaskTimes( "../out/", tasks );

  writeEVEvolution( "../out/gamma.dat", "../out/omega.dat", tasks );

  delete loadmodel;

  MPI_Finalize();
}

void
writeGroupTimes( const char* filename, ProcessGroupManagerContainer& pgroups )
{
  std::ofstream myfile( filename );

  // write first line with column descriptions
  myfile << "% \t total group time \t number of tasks" << std::endl;

  // write data
  for( size_t i = 0; i < pgroups.size(); ++i ){
    const TaskIDContainer tasks = pgroups[i].getTaskIDContainer();

    double totaltime( 0.0 );
    for( size_t j = 0; j < tasks.size(); ++j ){
      totaltime += tasks[i]->getTotalWallTime();
    }

    double ivtime( 0.0 );
    for( size_t j = 0; j < tasks.size(); ++j )
      ivtime += tasks[i]->getTimeIV()[0];

    myfile << i << " " << totaltime << " " << ivtime << " " << tasks.size()
    << std::endl;
  }

  myfile.close();
}

void
writeOverallTime( const char* filename, combigrid::real time, size_t ngroups,
  size_t ntasks )
{
  std::ofstream myfile( filename );

  // write first line with column descriptions
  myfile << "% number of groups \t overall time \t number of tasks"
  << std::endl;

  // write data
  myfile << ngroups << " " << time << " " << ntasks << std::endl;
}

void
writeSolutionTimes( const char* filename, GeneTaskContainer& tasks )
{
  std::ofstream myfile( filename );

  // write first line with column descriptions
  myfile << "% id \t level vector \t run time \t perf time" << std::endl;

  // write data
  for( size_t i = 0; i < tasks.size(); ++i ){
    LevelVector l = tasks[i].getLevelVector();
    combigrid::real time = tasks[i].getLastWallTime();
    int id = tasks[i].getID();
    combigrid::real time_perf = tasks[i].getTimePerf()[0];
    myfile << id << " " << l << " " << time << " " << time_perf << std::endl;
  }
}

void
writeTaskTimes( const char* outputdir, GeneTaskContainer& tasks )
{
  for( size_t i = 0; i < tasks.size(); ++i ){
    GeneTask& t = tasks[i];

    std::string filename = outputdir + std::string( "/task" )
    + boost::lexical_cast<std::string>( i ) + std::string( ".dat" );
    std::ofstream ofs( filename.c_str() );

    ofs << "%step" << " " << "time_autopar" << " " << "time_IV" << " "
    << "time_cp_mem" << " " << "time_cp" << " " << "time_gatherCP" << " "
    << "time_convertCPtoFG" << " " << "time_hierarchization" << " "
    << "time_addHFGtoSG" << " " << "time_reduce" << " "
    << "time_convertSGtoHFG" << " " << "time_dehierarchization" << " "
    << "time_convertFGtoCP" << " " << "time_scatterCP" << std::endl;

    for( size_t j = 0; j < t.getTimeIV().size(); ++j ){
      combigrid::real time_autopar( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimePerf();
        if( j < tmp.size() ){
          time_autopar = tmp[j];
        }
      }

      combigrid::real time_IV( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeIV();
        if( j < tmp.size() )
          time_IV = tmp[j];
      }

      combigrid::real time_cp_mem( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeCPMem();
        if( j < tmp.size() )
          time_cp_mem = tmp[j];
      }

      combigrid::real time_cp( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeCP();
        if( j < tmp.size() )
          time_cp = tmp[j];
      }

      combigrid::real time_gatherCP( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeGatherCP();
        if( j < tmp.size() )
          time_gatherCP = tmp[j];
      }

      combigrid::real time_convertCPtoFG( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeConvertCPtoFG();
        if( j < tmp.size() )
          time_convertCPtoFG = tmp[j];
      }

      combigrid::real time_hierarchization( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeHierarchization();
        if( j < tmp.size() )
          time_hierarchization = tmp[j];
      }

      combigrid::real time_addHFGtoSG( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeAddHFGtoSG();
        if( j < tmp.size() )
          time_addHFGtoSG = tmp[j];
      }

      combigrid::real time_reduce( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeReduce();
        if( j < tmp.size() )
          time_reduce = tmp[j];
      }

      combigrid::real time_convertSGtoHFG( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeConvertSGtoHFG();
        if( j < tmp.size() )
          time_convertSGtoHFG = tmp[j];
      }

      combigrid::real time_dehierarchization( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeDehierarchization();
        if( j < tmp.size() )
          time_dehierarchization = tmp[j];
      }

      combigrid::real time_convertFGtoCP( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeConvertFGtoCP();
        if( j < tmp.size() )
          time_convertFGtoCP = tmp[j];
      }

      combigrid::real time_scatterCP( -1.0 );
      {
        const std::vector<combigrid::real>& tmp = t.getTimeScatterCP();
        if( j < tmp.size() )
          time_scatterCP = tmp[j];
      }

      ofs << j << " " << time_autopar << " " << time_IV << " " << time_cp_mem
      << " " << time_cp << " " << time_gatherCP << " " << time_convertCPtoFG
      << " " << time_hierarchization << " " << time_addHFGtoSG << " "
      << time_reduce << " " << time_convertSGtoHFG << " "
      << time_dehierarchization << " " << time_convertFGtoCP << " "
      << time_scatterCP << std::endl;
    }
  }
}

void writeEVEvolution( const char* filename_gamma, const char* filename_omega,
 GeneTaskContainer& tasks ){
  if( tasks.size() == 0 )
    return;

  std::ofstream gammafile( filename_gamma );
  std::ofstream omegafile( filename_omega );

  // get times of first task
  std::vector<size_t>& steps = tasks[0].getLambdaTimeSteps();

  for( size_t i=0; i<steps.size(); ++i ){
    size_t timestep = steps[i];

    gammafile << timestep;
    omegafile << timestep;

    for( size_t j=0; j<tasks.size(); ++ j ){
      GeneTask& t = tasks[j];

      assert( t.getLambdaTimeSteps()[i] == timestep );
      gammafile << " " << t.getLambdaEvolution()[i].real();
      omegafile << " " << t.getLambdaEvolution()[i].imag();
    }

    gammafile << std::endl;
    omegafile << std::endl;
  }
}
