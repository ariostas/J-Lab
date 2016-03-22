for x in *.lvm;
do
   mv $x "${x:0: -4}.txt"
done
