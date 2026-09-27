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
    echo "Usage: $0 VERSION [--no-native] [--no-perf-native] [--no-system] [--editable]"
    echo "  --no-native   disable -march=native (default: enabled)"
    echo "  --no-perf-native use system packages for the performance shortlist"
    echo "  --no-system   build all dependencies from source"
    echo "  --editable    keep Sage in editable mode (default: disabled)"
    exit 1;
fi

version=$1; shift
USE_NATIVE=true
PERF_NATIVE=true
NO_SYSTEM=false
USE_EDITABLE=false
for arg in "$@"; do
    case "$arg" in
        --no-native)  USE_NATIVE=false ;;
        --no-perf-native) PERF_NATIVE=false ;;
        --no-system)  NO_SYSTEM=true ;;
        --editable)   USE_EDITABLE=true ;;
        *) echo "Unknown option: $arg"; exit 1 ;;
    esac
done

echo "Compiling version $version"
$USE_NATIVE && echo "  -march=native enabled"
$PERF_NATIVE && echo "  performance-critical packages built from source (default)"
! $PERF_NATIVE && echo "  using system packages for the performance shortlist"
$NO_SYSTEM && echo "  building all dependencies from source"
$USE_EDITABLE && echo "  editable Sage install enabled"
! $USE_EDITABLE && echo "  non-editable Sage install enabled"

echo "Sleeping 5s"
sleep 5s;

# run prerequisites_sage.sh first for system dependencies
archive="sage-${version}.tar.gz"
archive_url="https://mirrors.mit.edu/sage/src/${archive}"

if [ -f "$archive" ] && tar tzf "$archive" >/dev/null 2>&1; then
    echo "Using existing $archive"
else
    if [ -f "$archive" ]; then
        echo "Existing $archive is incomplete or invalid; attempting to resume download"
    fi
    wget -c "$archive_url"
    if ! tar tzf "$archive" >/dev/null 2>&1; then
        echo "Resumed download is still invalid; re-downloading $archive from scratch"
        rm -f "$archive"
        wget "$archive_url"
    fi
fi

if ! tar tzf "$archive" >/dev/null 2>&1; then
    echo "Downloaded archive $archive is invalid"
    exit 1
fi

tar xf "$archive"
cd sage-${version}

CONFIGURE_FLAGS=()
PERF_NATIVE_PACKAGES=(
    openblas
    gmp
    ntl
    fflas_ffpack
    linbox
    givaro
    m4ri
    m4rie
    ecm
    primecount
    pari
    flint
)

append_config_flag() {
    local flag=$1
    local existing
    for existing in "${CONFIGURE_FLAGS[@]}"; do
        if [ "$existing" = "$flag" ]; then
            return
        fi
    done
    CONFIGURE_FLAGS+=("$flag")
}

if ! $USE_EDITABLE; then
    if ./configure --help | grep -q -- '--disable-editable'; then
        append_config_flag --disable-editable
    else
        echo "  warning: this Sage version does not support --disable-editable; continuing without it"
    fi
fi

if $PERF_NATIVE && ! $NO_SYSTEM; then
    for pkg in "${PERF_NATIVE_PACKAGES[@]}"; do
        append_config_flag "--without-system-$pkg"
    done
fi

if $NO_SYSTEM; then
    while IFS= read -r flag; do
        append_config_flag "$flag"
    done < <(
        ls build/pkgs/*/spkg-configure.m4 | \
        sed 's|build/pkgs/\(.*\)/spkg-configure.m4|--without-system-\1|'
    )
fi

if $USE_NATIVE; then
    export CFLAGS="${CFLAGS} -march=native"
    export CXXFLAGS="${CXXFLAGS} -march=native"
    export FCFLAGS="${FCFLAGS} -march=native"
fi

./configure "${CONFIGURE_FLAGS[@]}"
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
./sage -pip install --no-build-isolation --upgrade git+https://github.com/edgarcosta/pyrforest.git
./sage -pip install --upgrade git+https://github.com/edgarcosta/pycontrolledreduction.git@master#egg=pycontrolledreduction
./sage -b
cd ..
chmod a+rX -R sage-${version}
echo "If you want this to be the new version to be used, don't forget to do:"
echo "$ rm ~/sage-root && ln -s sage-${version} ~/sage-root"
set +e
