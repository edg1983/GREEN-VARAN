read -r -d '' HELP_MESSAGE << EOM
### PREPARE SCORES UTILITY ###
Author: Edoardo Giacopuzzi

Utility to download and prepare non-coding prediction scores
used by GREEN-VARAN for annotation

Usage:
    build_score -s score -b [GRCh37 | GRCh38] [-o output_folder] [-h]

Output folder default to /scores within the current directory 
Allowed score values:
CADD, DANN, FIRE, LinSight, NCBoost, ReMM
- all: download and process all scores
- essential: download ReMM, LinSight, NCBoost
EOM

ALLSCORES=(CADD DANN FIRE LinSight NCBoost ReMM)
ESSENTIALS=(LinSight NCBoost ReMM)
LIFTOVER="liftOver"
CHAIN="###"

declare -A REPOS
REPOS=(
    ["CADD_GRCh38"]="https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz" \
    ["CADD_GRCh37"]="https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz" \
    ["DANN_GRCh37"]="https://cbcl.ics.uci.edu/public_data/DANN/data/DANN_whole_genome_SNVs.tsv.bgz" \
    ["FIRE_GRCh37_chr1"]="https://www.dropbox.com/s/fpxxbzrbi460te7/FIRE_chr1.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr2"]="https://www.dropbox.com/s/6fb97os3jjik91z/FIRE_chr2.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr3"]="https://www.dropbox.com/s/prljbprygiis08o/FIRE_chr3.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr4"]="https://www.dropbox.com/s/p03h2zqc8ozl98y/FIRE_chr4.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr5"]="https://www.dropbox.com/s/vbhhc35ae64i5ei/FIRE_chr5.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr6"]="https://www.dropbox.com/s/ro2ay1rp6o4elwt/FIRE_chr6.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr7"]="https://www.dropbox.com/s/dxi23frpt4iwcmh/FIRE_chr7.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr8"]="https://www.dropbox.com/s/yzwjyf3pmbp6zct/FIRE_chr8.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr9"]="https://www.dropbox.com/s/na8zrxajovhmxvn/FIRE_chr9.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr10"]="https://www.dropbox.com/s/x21ifcd6vohud1a/FIRE_chr10.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr11"]="https://www.dropbox.com/s/bl791xo01s03f3n/FIRE_chr11.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr12"]="https://www.dropbox.com/s/tqd1pcv88eh2dr8/FIRE_chr12.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr13"]="https://www.dropbox.com/s/0rzzfw90974jesk/FIRE_chr13.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr14"]="https://www.dropbox.com/s/zox5qin56qqjzi0/FIRE_chr14.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr15"]="https://www.dropbox.com/s/g8k0e37vjd3q9j0/FIRE_chr15.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr16"]="https://www.dropbox.com/s/w0wvnovn2ramr46/FIRE_chr16.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr17"]="https://www.dropbox.com/s/yex2381mw0c34c3/FIRE_chr17.txt.gz?dl=1" \
    ["FIRE_GRCh37_chr18"]="http://web.stanford.edu/~nilah/fire/FIRE_chr18.txt.gz" \
    ["FIRE_GRCh37_chr19"]="http://web.stanford.edu/~nilah/fire/FIRE_chr19.txt.gz" \
    ["FIRE_GRCh37_chr20"]="http://web.stanford.edu/~nilah/fire/FIRE_chr20.txt.gz" \
    ["FIRE_GRCh37_chr21"]="http://web.stanford.edu/~nilah/fire/FIRE_chr21.txt.gz" \
    ["FIRE_GRCh37_chr22"]="http://web.stanford.edu/~nilah/fire/FIRE_chr22.txt.gz" \
    ["FIRE_GRCh37_chrX"]="https://www.dropbox.com/s/12pylu6tixvcw6p/FIRE_chrX.txt.gz?dl=1" \
    ["FIRE_GRCh37_chrY"]="https://www.dropbox.com/s/x1he9l1byw50wzk/FIRE_chrY.txt.gz?dl=1"
)

if [ $# == 0 ]
then
    echo -e "Missing mandatory arguments\n"
    echo "$HELP_MESSAGE"
    exit 1
fi

if [ -x "$DIR/liftOver" ]
then
    echo "liftOver file not found / executable in $DIR"
fi

if [ ! -f "$DIR/hg19toGRCh38.chain.gz" ] 
then
    echo "liftOver chain file hg19toGRCh38.chain.gz not found in $DIR"
fi

while test $# -gt 0; do
    case "$1" in
        -h|--help)
            echo "$HELP_MESSAGE"
            exit 0
        ;;
        -s|--score)
            shift
            if test $# -gt 0; then
                export SCORE=$1
            else
                echo "No score name specified. Use -h for help"
                exit 1
            fi
            shift
        ;;
        -b|--build)
            shift
            if test $# -gt 0; then
                if [ $1 == "GRCh37" ] || [ $1 == "GRCh38" ]; then
                    export BUILD=$1
                else
                    echo "One of GRCh37 or GRCh38 must be specified with -b/--build"
                    exit 1
                fi
            else
                echo "No genome build specified. Use -h for help"
                exit 1
            fi
            shift
        ;;        
        -o|--out)
            shift
            if test $# -gt 0; then
                export OUTDIR=$1
            else
                export OUTDIR="scores"
            fi
            shift
        ;;
        -l|--liftover)
            shift
            if test $# -gt 0; then
                export LIFTOVER=$1
            else
                echo "Please provide location of liftOver tool"
            fi
            shift
        ;;
        -c|--chain)
            shift
            if test $# -gt 0; then
                export CHAIN=$1
            else
                echo "Please provide location of hg19Tohg38 chain file"
            fi
            shift
        ;;
        *)
            break
        ;;
    esac
done

liftover_bed () {
    local in=$1
    local out=$2
    $LIFTOVER $in $CHAIN TMP_GRCh38 TMP_dropped
    echo "#HEADER" > $out
    cat TMP_GRCh38 | tr "@" "\t" | cut -f1,3- >> $out
}

make_bed () {
    local in=$1
    local out=$2
    local chr_prefix=$3
    local format=$4
    local sep=$5
    local startline=$6
    if [ format == "tsv" ]; then
        paste -d "\t" \
        <(zcat $in | tail -n+$startline | awk -F"$sep" -v chr="$chr_prefix" '{OFS="\t"}; {print chr$1, $2-1, $2}') \
        <(zcat $in | tail -n+$startline | cut -f3- | tr "$sep" "@") \
        > $out
    fi
    if [ format == "bed" ]; then
        paste -d "\t" \
        <(zcat $in | tail -n+$startline | awk -F"$sep" -v chr="$chr_prefix" '{OFS="\t"}; {print $1, $2, $3}') \
        <(zcat $in | tail -n+$startline | cut -f4- | tr "$sep" "@") \
        > $out
    fi
}

download () {
    local myscore=$1
    echo "Downloading $myscore annotations"
    case $myscore in
        CADD)
            wget -nv "${REPOS[${myscore}_GRCh37]}" -O TMP_GRCh37_${myscore}.tsv.gz
            wget -nv "${REPOS[${myscore}_GRCh37]}" -O TMP_GRCh38_${myscore}.tsv.gz
        ;;
        DANN)
            wget -nv "${REPOS[${myscore}_GRCh37]}" -O TMP_GRCh37_${myscore}.tsv.gz
            make_bed(TMP_GRCh37_${myscore}.tsv.gz, TMP_GRCh37_${myscore}.bed, "chr", "tsv", "\t", 2) 
            liftover_bed(TMP_GRCh37_${myscore}.bed, TMP_GRCh38_${myscore}.tsv.gz)
        ;;
        FIRE)
            for c in {1..22} X Y; do
                wget -nv "${REPOS[${myscore}_GRCh37_chr${c}]}" -O TMP_GRCh37_${myscore}_chr${c}.tsv.gz
            done
            head -1 ${REPOS[${myscore}_GRCh37_chr1]} > TMP_GRCh37_${myscore}.tsv
            for c in {1..22} X Y; do
                zcat TMP_GRCh37_${myscore}_chr${c}.tsv.gz | tail -n+2 >> TMP_GRCh37_${myscore}.tsv
            done
            bgzip TMP_GRCh37_${myscore}.tsv
            make_bed(TMP_GRCh37_${myscore}.tsv.gz, TMP_GRCh37_${myscore}.bed, "chr", "tsv", "\t", 2)
            liftover_bed(TMP_GRCh37_${myscore}.bed, TMP_GRCh38_${myscore}.tsv.gz)
        ;;
    esac
}

compress () {
    local myscore=$1
    local s=$2
    local b=$3
    local e=$4
    bgzip GRCh37_${myscore}.tsv
    bgzip GRCh38_${myscore}.tsv
    tabix GRCh37_${myscore}.tsv -s $s -b $b -e $e -m 12 --csi -f
    tabix GRCh38_${myscore}.tsv -s $s -b $b -e $e -m 12 --csi -f
}

preprocess () {
    local myscore=$1
    case $myscore in
        CADD)
            echo -e "#CHROM\tPOS\tREF\tALT\tSCORE" > GRCh37_${myscore}.tsv
            echo -e "#CHROM\tPOS\tREF\tALT\tSCORE" > GRCh38_${myscore}.tsv
            zcat TMP_GRCh37_${myscore}.tsv.gz | grep -v "#" | cut -f1-4,6 >> GRCh37_${myscore}.tsv & \
            zcat TMP_GRCh38_${myscore}.tsv.gz | grep -v "#" | cut -f1-4,6 >> GRCh38_${myscore}.tsv
        ;;
        DANN)
            echo -e "#CHROM\tPOS\tREF\tALT\tSCORE" > GRCh37_${myscore}.tsv
            echo -e "#CHROM\tPOS\tREF\tALT\tSCORE" > GRCh38_${myscore}.tsv
            zcat TMP_GRCh37_${myscore}.tsv.gz | grep -v "#" | awk -F"\t" '{OFS="\t"}; NR > 1 {$5=substr($5,1,5)}; {print;}' >> GRCh37_${myscore}.tsv & \
            zcat TMP_GRCh37_${myscore}.tsv.gz | grep -v "#" | awk -F"\t" '{OFS="\t"}; NR > 1 {$5=substr($5,1,5)}; {print;}' >> GRCh37_${myscore}.tsv
        ;;
        FIRE)
            echo -e "#CHROM\tPOS\tREF\tALT\tSCORE" > GRCh37_${myscore}.tsv
            echo -e "#CHROM\tPOS\tREF\tALT\tSCORE" > GRCh38_${myscore}.tsv
            zcat TMP_GRCh37_${myscore}.tsv.gz | tail -n+2 >> GRCh37_${myscore}.tsv & \
            zcat TMP_GRCh38_${myscore}.tsv.gz | tail -n+2 >> GRCh38_${myscore}.tsv   
        ;;
    esac
}

#Check commands and files
if ! command -v bgzip >/dev/null 2>&1 || ! command -v tabix >/dev/null 2>&1 || ! command -v wget >/dev/null 2>&1 || ! command -v $LIFTOVER >/dev/null 2>&1; then
    echo "tabix, bgzip, wget or liftover command not executable!!"
    exit 1
fi

if [ ! -f $CHAIN ]; then
    echo "Chain file $CHAIN not found!!"
    exit 1
fi

#Start score download and process
echo "Output folder: $OUTDIR"
mkdir -p $OUTDIR
cd $OUTDIR

case $SCORE in 
    all)
        echo "All scores selected... This can take a long time..."
        for myscore in "${ALLSCORES[@]}"; do 
            echo "# Start $myscore"
            filename="${BUILD}_${myscore}"
            raw_file=$(download($myscore, "${REPOS[${myscore}_$BUILD]}", "$filename"))
            processed_file=$(preprocess($myscore, "$raw_file", "$filename"))
            compress("$processed_file.tsv")
            rm $raw_file
        done
    ;;
    essential) 
        echo "Essential scores selected: LinSight, NCBoost, ReMM"
        for myscore in "${ESSENTIALS[@]}"; do 
            echo "# Start $myscore"
            filename="${BUILD}_${myscore}"
            raw_file=$(download($myscore, "${REPOS[${myscore}_$BUILD]}", "$filename"))
            processed_file=$(preprocess($myscore, "$raw_file", "$filename"))
            compress("$processed_file.tsv")
            rm $raw_file
        done
    ;;
    test)
        echo "test complete"
    ;;
    *)
        download($SCORE)
        preprocess($SCORE)
        compress($SCORE)
        rm $raw_file
    ;;
esac

rm TMP_*
cd -
echo "Finished processing!"