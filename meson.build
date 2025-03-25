# This file is part of QCxMS2.

project(
  'qcxms2',
  'fortran',
  version: '1.0.0',
  license: 'GPL-3.0-or-later',
  meson_version: '>=0.55',
  default_options: [
    'buildtype=release',
    'default_library=static',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
exe_deps = []
subdir('config')

srcs = []
prog = []
# Collect sources
subdir('src')

# Create library target
qcxms_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: exe_deps,
  #  include_directories: include_directories('subprojects'),
)

# Export as dependency
qcxms_inc = qcxms_lib.private_dir_include()
qcxms_dep = declare_dependency(
  link_with: qcxms_lib,
  include_directories: qcxms_inc,
  dependencies: exe_deps,
)


commit = get_option('build_name')
git = find_program('git', required: false)
if git.found()
  git_commit = run_command(git, 'show', '-s', '--format=%h')
  if git_commit.returncode() == 0
    commit = git_commit.stdout().strip()
  endif
endif
# We rely on the existence of Python since meson will need a Python interpreter to run,
# this way we have a cross-platform way to obtain relevant data for our build.
pymod = import('python')
python = pymod.find_installation('python3', required: false)
if not python.found()
  python = find_program('python3', 'python')
endif

# build a configuration data containing all the important data to propagate
# it to the automatically generated files
config = configuration_data({
  'name': meson.project_name(),
  'description': 'QM-based Mass Spectra computation package',
  'version': meson.project_version(),
  'commit': commit,
  'date': run_command(python, '-c', 'import datetime; print(datetime.date.today().strftime("%Y-%m-%d"))').stdout().strip(),
  'author': run_command(python, '-c', 'import getpass; print(getpass.getuser())').stdout().strip(),
  'origin': run_command(python, '-c', 'import socket; print(socket.gethostname())').stdout().strip(),
  'prefix': get_option('prefix'),
  'mandir': get_option('mandir'),
  'bindir': get_option('bindir'),
  'libdir': get_option('libdir'),
  'includedir': get_option('includedir'),
  'datadir': get_option('datadir')/meson.project_name(),
})

configure_file(
  input: files('assets/templates/version.f90'),
  output: 'qcxms2_version.fh',
  configuration : config,
)

qcxms2_header = files('include/qcxms2.h')

## ========================================== ##
## INSTALL
## ========================================== ##
if install
  install_headers(qcxms2_header)
endif


# Create executable target
qcxms2_exe = executable(
  meson.project_name(),
  sources: prog,
  dependencies: qcxms_dep,
  install: install,
)
