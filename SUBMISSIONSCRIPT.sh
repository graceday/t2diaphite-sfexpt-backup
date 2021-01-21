#! /bin/bash
#$ -S /bin/bash
#$ -l h=!comp106&!comp143&!comp144&!comp145&!comp146&!comp147&!comp148&!comp149&!comp150&!comp151
#$ -l s_rt=150:00:00
#$ -cwd
#$ -M grace.day@new.ox.ac.uk
#$ -m eas

# Use this variable to track if copying files has worked successfully.
# Leave the files on the node for manual intervention if it has not.
SUCCESSFLAG="true"
NODENAME=$(uname -n)
STARTTIME=$(date '+%Y-%m-%d-%H-%M')
#If we are running this script locally, then the JOB_ID isn't set.
# Set it here as a temporary date.
if [ -z "$JOB_ID" ]; then
    JOB_ID=$(date '+%Y-%m-%d-%H-%M')
fi
echo "Started at $STARTTIME on node $NODENAME with JOBID $JOB_ID."

# FROMDIR is the directory your input files are in, and TODIR is where they're put afterwards
FROMDIR=$(pwd)
TODIR=$(pwd)


# ALTERNATIVE TODIR, date style:
#dt=$(date '+%Y-%m-%d-%H-%M')
#TODIR="`pwd/Job-$dt"

# INFILES are the the files we want to copy across as inputs.
declare -a INFILES=(*.data *.inpt)
# OUTFILES are the files we want to copy back home.
declare -a OUTFILES=(*.lammpstrj log.* *.dat)
# WORKDIR is a scratch directory, cleared up at the end since $SCRATCH
# isn't always defined, and EXEC is the name of the executable
WORKDIR="/tmp/$USER/$JOB_ID"
EXEC="lmp_serial"
ARGS="-in *.inpt"
echo "Executing $EXEC in work folder $WORKDIR. Copying back to $TODIR."

export WORKDIR
mkdir -p "$WORKDIR"
mkdir -p "$WORKDIR/output_files/"
mkdir -p "$TODIR"

cd "$FROMDIR"
for FILE in ${INFILES[@]}; do
    cp "$FILE" "$WORKDIR"
done

# Move into the work directory and execute the program.
cd "$WORKDIR"
"$EXEC" $ARGS

# OUTFILES get copied across to TODIR when you're done
# ZIPFILES get compressed before sending over.
if [ ${#OUTFILES[@]} -eq 0 ]; then
    echo "Could not find any output files."
    SUCCESSFLAG="false"
fi


declare -a ZIPFILES=(*.lammpstrj)
for FILE in ${ZIPFILES[@]}; do
    if [[ -f "$FILE" ]]; then
        echo "Gzipping $FILE"
        gzip "$FILE"
        OUTFILES+=($FILE.gz)
    else
        echo "Failed to gzip $FILE -- it appears to have gone missing." 1>&2
        SUCCESSFLAG="false"
    fi
done

for FILE in ${OUTFILES[@]}; do
    if [[ -f "$FILE" ]]; then
        cp "$FILE" "$TODIR"
    else
        if [[ " ${ZIPFILES[@]} " =~ " ${FILE} " ]]; then
            # Phew, this file has gone missing thanks to gzip.
            echo "$FILE has been zipped."  1>&2
        else
            echo "Failed to copy $FILE -- it appears to have gone missing." 1>&2
            SUCCESSFLAG="false"
        fi     
    fi
done

cd "$TODIR"

declare -a FOUNDFILES=("$TODIR"/*)
if [ ${#OUTFILES[@]} -eq 0 ]; then
    echo "Could not find any output files transferred over."
    SUCCESSFLAG="false"
fi

# Now unzip all of the files again.
gunzip "$TODIR/"*.gz
# Tidy up after ourselves in the scratch
ENDTIME=$(date '+%Y-%m-%d-%H-%M')
echo "Job terminated successfully at $ENDTIME on node $NODENAME with JOBID $JOB_ID."  1>&2
if [[ SUCCESSFLAG = "false" ]]; then
    
    echo "Job has not tidied up after itself. Please check $NODENAME:$WORKDIR and retrieve the files manually."  1>&2
else
    rm -Rf "$WORKDIR"
fi
