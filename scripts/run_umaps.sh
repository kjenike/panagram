

anchor=$1
fasta=$2

window=100000
eps=2
n_neighbor=3
md=0
skips=100


#for i in {1..38};
while read col1 col2 colmore
do
	echo "$col1"

	python plot_umap.py $anchor "$col1" $window "$eps" "$n_neighbor" "$md" "$skips" > anchor/$anchor/umap_clusters_"$col1"_"$window"_"$eps"_"$n_neighbor"_"$md"_"$skips".txt

done < $fasta.fai


