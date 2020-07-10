counter=1
while [ $counter -le 3 ]
do

$1 --bfile $2 --make-grm-part 3 $counter --out $3	

((counter++))
done