#!/bin/bash
_user="$(id -u -n)"
_uid="$(id -u)"
echo "User name : $_user"
echo "User name ID (UID) : $_uid"
j=$(nproc 2>/dev/null || echo 4)
echo "Using $j threads to compile"
set -e
#uid = 1200
if [ $_uid -ne 1200 ]
then
    echo "Only the user sage with UID = 1200 should be installing sage"; exit 1;
fi

if [ -z "$1" ]
then
    echo "Usage: $0 VERSION [--native] [--no-system]"
    echo "  --native      compile with -march=native (AVX-512, etc.)"
    echo "  --no-system   build all dependencies from source"
    exit 1;
fi

version=$1; shift
USE_NATIVE=false
NO_SYSTEM=false
for arg in "$@"; do
    case "$arg" in
        --native)    USE_NATIVE=true ;;
        --no-system) NO_SYSTEM=true ;;
        *) echo "Unknown option: $arg"; exit 1 ;;
    esac
done

echo "Compiling version $version"
$USE_NATIVE && echo "  -march=native enabled"
$NO_SYSTEM && echo "  building all dependencies from source"

echo "Sleeping 5s"
sleep 5s;

# run prerequisites_sage.sh first for system dependencies
wget https://mirrors.mit.edu/sage/src/sage-${version}.tar.gz -O sage-${version}.tar.gz
tar xf sage-${version}.tar.gz
cd sage-${version}

CONFIGURE_FLAGS="--without-system-python3"

if $NO_SYSTEM; then
    CONFIGURE_FLAGS=$(ls build/pkgs/*/spkg-configure.m4 | \
        sed 's|build/pkgs/\(.*\)/spkg-configure.m4|--without-system-\1|' | tr '\n' ' ')
fi

if $USE_NATIVE; then
    export CFLAGS="${CFLAGS} -march=native"
    export CXXFLAGS="${CXXFLAGS} -march=native"
    export FCFLAGS="${FCFLAGS} -march=native"
fi

./configure $CONFIGURE_FLAGS
MAKE="make -j${j}" make
./sage -i gap_packages
./sage -b
wget https://raw.githubusercontent.com/LMFDB/lmfdb/main/requirements.txt
./sage -pip install -r requirements.txt
wget https://raw.githubusercontent.com/roed314/seminars/master/requirements.txt -O semrequirements.txt
./sage -pip install -r semrequirements.txt --upgrade
wget https://raw.githubusercontent.com/AndrewVSutherland/psetpartners/master/requirements.txt -O psetrequirements.txt
./sage -pip install -r psetrequirements.txt --upgrade
./sage -pip install bcrypt
./sage -pip install gunicorn pyflakes
./sage -pip install greenlet eventlet gevent
./sage -b
cd ..
chmod a+rX -R sage-${version}
echo "If you want this to be the new version to be used, don't forget to do:"
echo "$ rm ~/sage-root && ln -s sage-${version} ~/sage-root"
set +e
