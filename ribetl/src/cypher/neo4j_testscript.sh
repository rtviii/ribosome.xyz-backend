#!/usr/bin/bash
set -e
# set -u
set -o pipefail

# For remote access:
# echo "match (n) return count(n)" | cypher-shell -a neo4j+s://ribosome.xyz:7687 --format plain -u 'rt' -p 'rrr' --database 'riboauth'


usage (){
      echo "script usage: $0 [-f] [path to strucutres' json profile] [-r] [neo4j remote server i.e. localhost:7678]" >&2
      exit 1
}


DEFAULT_NEO4J_ADDRESS="127.0.0.1:7687"
while getopts 'd:f:a:' OPTION; do
  case "$OPTION" in
    d)
      DATABASE_NAME="$OPTARG"
      ;;
    f)
      FILEPATH="$OPTARG"
      if [ -f $FILEPATH ];
      then
        file=$FILEPATH
        extension=${file##*.}
        if [ $extension != "json" ];
        then
          echo "The profile file must be a .json. Exiting."
          exit $((2))
        fi
        structid=${file::4}
        structid=${structid^^}
              echo "Processing $structid"
      else
        echo "$FILEPATH is not an acceptable file"
        exit $((1))
      fi
      ;;
    a)
      NEO4J_ADDRESS="$OPTARG"
      ;;
    ?)
    usage
      ;;
  esac
done

if [ -z $NEO4J_ADDRESS ];
then
  NEO4J_ADDRESS="$DEFAULT_NEO4J_ADDRESS"
fi


if [ -z $FILEPATH] || [ -z $DATABASE_NAME];
then
  usage
fi

CYPHER="match (n:RibosomeStructure) return count(n);"
echo $CYPHER | cypher-shell -a "$NEO4J_ADDRESS" --format plain -u 'rt' -p 'rrr' --database $DATABASE_NAME


# # echo "X into shell1"
# # echo "X into shell2" in zsh