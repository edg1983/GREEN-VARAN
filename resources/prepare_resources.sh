read -r -d '' HELP_MESSAGE << EOM
### PREPARE SCORES UTILITY ###
Author: Edoardo Giacopuzzi

Utility to download and prepare resources used by GREEN-VARAN tools.
The utility will download the required annotation and eventually
create expected sub-folder.

Usage:
    prepare_resources.sh
        [-r,--resource RESOURCES_FILE]  GREEN-VARAN_resources.txt (default)
        [-n,--name RESOURCE_NAME]       all = download all (default)
                                        scores = download all scores
                                        bed_files = download all bed_files
                                        AF = download AF files
                                        SV_annotation = download all SV annotations
                                        specific_name (see --list)
        [-o,--out OUTPUT_DIR]           current directory (default)
        [-l,--list]                     list available resources
        [-h,--help]                     get help message
EOM

function downloadRepoGroup {
    repo_class = $1
    n_repos=$(cat $RESOURCE_FILE | awk -F"\t" -v name="$repo_class" '{$2 == name}' | wc -l)
        echo "$n_repos scores found in resource file"
        while read -r repo_id repo_dir repo_http; do
            mdkir -p $repo_dir
            cd $repo_dir
            echo "Downloading $repo_id in $OUTDIR/$repo_dir"
            wget --backups=1 -nv $repo_http
            cd ..    
        done < $(awk -F"\t" -v name="$repo_name" '{$2 == name}')   
}

if [ $# == 0 ]
then
    echo -e "Missing mandatory arguments\n"
    echo "$HELP_MESSAGE"
    exit 1
fi

while test $# -gt 0; do
    case "$1" in
        -h|--help)
            echo "$HELP_MESSAGE"
            exit 0
        ;;
        -r|--resource)
            shift
            if test $# -gt 0; then
                export RESOURCE_FILE=$1
            else
                echo "No resource file specified (-r)"
                echo "Using default: GREEN-VARAN_resources.txt"
                export RESOURCE_FILE="GREEN-VARAN_resources.txt"
            fi
            shift
        -n|--name)
            shift
            if test $# -gt 0; then
                    export REPO_NAME=$1
            else
                echo "No resource name specified (-n)"
                echo "Downloading all resources"
                export REPO_NAME="all"
            fi
            shift
        ;;        
        -o|--out)
            shift
            if test $# -gt 0; then
                export OUTDIR=$1
            else
                export OUTDIR="./"
            fi
            shift
        ;;
        -l|--list)
            export LIST_RES=1
            shift
        ;;
        *)
            break
        ;;
    esac
done

#Check resource file
if [ ! -f $RESOURCE_FILE ]; then
    echo "FATAL! - Resource file $RESOURCE_FILE do not exists!"
    exit 1
fi

#List resources if -l
if [ $LIST_RES == 1 ]; then
    echo "### LIST OF AVAILABLE RESOURCES ###"
    cut -f1 $RESOURCE_FILE
else
    echo "## OPTIONS AS INTERPRETED ##"
    echo -e "\tresource file: $RESOURCE_FILE
    output dir: $OUTDIR
    resource name: $REPO_NAME
    "
fi

#Check commands
if ! command -v gzip >/dev/null 2>&1 || ! command -v tar >/dev/null 2>&1 || ! command -v wget >/dev/null 2>&1; then
    echo "FATAL! - gzip, wget or tar command not executable!!"
    exit 1
fi

#Start download and process
current_dir=$(pwd)
mkdir -p $OUTDIR
cd $OUTDIR

case "$REPO_NAME" in
    all)
        n_repos=$(cat $RESOURCE_FILE | wc -l)
        echo "$n_repos sources found in resource file"
        while read -r repo_id repo_dir repo_http; do
            mdkir -p $repo_dir
            cd $repo_dir
            echo "Downloading $repo_id in $OUTDIR/$repo_dir"
            wget --backups=1 -nv $repo_http
            cd ..    
        done < $RESOURCE_FILE
    ;;
    scores|bed_files|SV_annotation|AF)
        downloadRepoGroup($REPO_NAME)
    ;;
    *)
        repo_id=$(grep $REPO_NAME $RESOURCE_FILE | cut -f1)
        repo_dir=$(grep $REPO_NAME $RESOURCE_FILE | cut -f2)
        repo_http=$(grep $REPO_NAME $RESOURCE_FILE | cut -f3)
        if [ "$repo_id" == ""]; then
            echo "FATAL! - Resource name $REPO_NAME not found in resource file. Use -l to list available resources"
            exit 1
        fi
        mdkir -p $repo_dir
        cd $repo_dir
        echo "Downloading $repo_id in $OUTDIR/$repo_dir"
        wget --backups=1 -nv $repo_http  
    ;;
esac

cd $current_dir
echo "All done!"