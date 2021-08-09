BASEDIR=$(pwd)
cd $BASEDIR

for dir in `ls -d */`
do 
    Kval=$(echo $dir | awk -F - '{print $NF}' | awk -F DPGP '{print $1}')
    for i in `ls $dir`:
    do
        #mv "$dir$i" "$Kval$i"
        file=$(echo $i | awk -F : '{print$1}')
        mv "$dir$file" "$Kval$file"
    done
done