#!/usr/bin/bash
set -e
set -u
set -o pipefail

# For remote access:
# echo "match (n) return count(n)" | cypher-shell -a neo4j+s://ribosome.xyz:7687 --format plain -u 'rt' -p 'rrr' --database 'riboauth'




# filepath=$1
# while getopts f:r: flag



# if [ -f $filepath ];
# then
# 	file=$(basename -- $filepath)
# 	extension=${file##*.}
# 	if [ $extension != "json" ];
# 	then
# 		echo "The profile file must be a .json. Exiting."
# 		exit $((2))
# 	fi
# 	structid=${file::4}
# 	structid=${structid^^}
#         echo "Processing $structid"
# else
# 	echo "$filepath is not an acceptable file"
# 	exit $((1))
# fi


unset -v remote
unset -v filepath

while getopts 'f:r:' OPTION; do
  case "$OPTION" in
    f)
      filepath="$OPTARG"
      echo "received $filepath"
      ;;
    r)
      remote="$OPTARG"
      echo "received $remote"
      ;;
    ?)
      echo "script usage: $0 [-f] [path to strucutres' json profile] [-r] [neo4j remote server i.e. localhost:7678]" >&2
      exit 1
      ;;
  esac
done

shift "$(($OPTIND -1))"

if [ -z "$filepath" ] ; then
        echo 'Provide path to .json profile of the structure. Exiting' >&2
        exit 1
fi

# echo "X into shell1"
# echo "X into shell2"