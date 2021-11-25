#Compile a static build with no dependencies on hts lib
#See: https://github.com/brentp/hts-nim#static-builds
#local '/well/gel/HICF2/software/singularity/musl-hts-nim.sif' 
pkg_name="greenvaran"

singularity exec \
--bind $PWD:$PWD \
--bind $PWD:/load/ \
docker://brentp/musl-hts-nim:latest \
/usr/local/bin/nsb \
-s $PWD/src/${pkg_name}.nim \
--nimble-file $PWD/${pkg_name}.nimble -- -d:danger -d:release
