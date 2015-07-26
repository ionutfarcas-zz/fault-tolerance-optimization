import os

inc = 'inc/'
build = 'build/'
src = 'src/'
lib = 'lib/'

lib_name = 'lp_opt'
target = 'opt_fault_tol'

if not os.path.exists(build):
	os.makedirs(build)

if not os.path.exists(lib):
	os.makedirs(lib)

env = Environment(CC = 'g++', ENV = os.environ)

env.Append(CPPPATH = ['inc/'])
env.Append(LIBPATH = ['build/'])
env.Append(SRCPATH = ['src/'])

env.Append(LIBS = ['-lm', '-lglpk'])
env.Append(CCFLAGS= ['-Wall', '-Werror', '-O2', '-std=c++11'])

env.VariantDir(variant_dir = build, src_dir = src, duplicate = 0)

src_files = Glob(build + '*.cpp')
env.Program(target, src_files)

env.SharedLibrary(target = lib + lib_name, source=src_files)
env.StaticLibrary(target = lib + lib_name, source=src_files)

