for data in 0, 1
do

for alg in 0, 1, 2, 3
do

#2, 5, 10, 12, 15
for eplison in 10
do

#1, 2, 3, 4, 5
for lambda in 1
do

for hour in 11
do

#20, 30, 40, 50, 60
for spread in 20, 30, 40, 50, 60
do


for monitor in 3600
do

#9, 12, 15, 18, 21
for duration in 15
do

for prune in 0
do

./tbip ${data} ${alg} ${eplison} ${lambda} ${hour} ${spread} ${monitor} ${duration} ${prune}

done
done
done
done
done
done
