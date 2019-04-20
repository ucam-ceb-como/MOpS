#!/bin/bash
#Prints out PAH structures.
#
#This script converts the dot files which are printed out using the --ppah flag
#into jpeg images. Note that you have to create the "JPG files" folder in the
#working directory.
#
#Written by Edward Yapp ekyy2@cam.ac.uk

sourcedots=KMC_DEBUG
resultjpgs=JPG\ files
for dotnames in `ls "$sourcedots"`
do
jpgname=${dotnames/.dot/.jpg}
dotfiles="$sourcedots"/"$dotnames"
jpgfiles="$resultjpgs"/"$jpgname"
neato -Tjpg "$dotfiles" -o "$jpgfiles"
echo "$dotnames"\ graph\ drawing\ Done!!
done

sourcedots=KMC_DEBUG/1
resultjpgs=JPG\ files/1
for dotnames in `ls "$sourcedots"`
do
jpgname=${dotnames/.dot/.jpg}
dotfiles="$sourcedots"/"$dotnames"
jpgfiles="$resultjpgs"/"$jpgname"
neato -Tjpg "$dotfiles" -o "$jpgfiles"
echo "$dotnames"\ graph\ drawing\ Done!!
done

#sourcedots=KMC_DEBUG/674
#resultjpgs=JPG\ files/674
#for dotnames in `ls "$sourcedots"`
#do
#jpgname=${dotnames/.dot/.jpg}
#dotfiles="$sourcedots"/"$dotnames"
#jpgfiles="$resultjpgs"/"$jpgname"
#neato -Tjpg "$dotfiles" -o "$jpgfiles"
#echo "$dotnames"\ graph\ drawing\ Done!!
#done

#sourcedots=KMC_DEBUG/784
#resultjpgs=JPG\ files/784
#for dotnames in `ls "$sourcedots"`
#do
#jpgname=${dotnames/.dot/.jpg}
#dotfiles="$sourcedots"/"$dotnames"
#jpgfiles="$resultjpgs"/"$jpgname"
#neato -Tjpg "$dotfiles" -o "$jpgfiles"
#echo "$dotnames"\ graph\ drawing\ Done!!
#done