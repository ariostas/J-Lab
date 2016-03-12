for x in *.Spe;
do
   sed -n 13,2060p $x > "${x:0: -4}.txt"
   rm $x
done
