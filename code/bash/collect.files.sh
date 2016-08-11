# Collecting files of interest

# This code is to rename and move files from our mothur
# folder to our dataEdited folder for subsequent analysis
# by other programs

cp data/WAVES/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared data/WAVES/dataEdited/WAVES.final.shared

cp data/WAVES/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list data/WAVES/dataEdited/WAVES.final.list

cp data/WAVES/mothur/data/WAVES/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy data/WAVES/dataEdited/WAVES.final.taxonomy

cp data/WAVES/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table data/WAVES/dataEdited/WAVES.final.count_table

cp data/WAVES/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.subsample.shared data/WAVES/dataEdited/WAVES.final.subsample.shared


# Then pull it off the Volume and into my computer
# Do this command in a terminal not connected to the instance.
# scp -i ~/Documents/gradSchool/research/edamameShellLesson.pem ubuntu@ec2-54-243-23-15.compute-1.amazonaws.com:/home/ubuntu/data/WAVES/dataEdited/* ~/Documents/gradSchool/research/mercury/WAVES/dataEdited/
