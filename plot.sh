for i in {0..9}
do
   cat ep.gp | sed 's/XXX/'$i'/g' > plotthis.gp
   gnuplot plotthis.gp
done
