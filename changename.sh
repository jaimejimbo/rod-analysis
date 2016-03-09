#!/bin/bash

#este directorio debe ser el lugar donde estén las fotos que quieres cambiar el nombre
FILES=*.JPG
i=1039
mkdir carpetaNueva

for f in $FILES
do
  echo "Procesando: $f"
  #no estoy seguro si tienes que poner ".png" o el tipo de foto que tengas al final del $i. Si no funciona, prueba cambiar por: "0_"$i".png" 
  if (($i<10)); then
    text="IMG_000"$i
  elif ((i<100)); then
    text="IMG_00"$i
  elif (($i<1000)); then
    text="IMG_0"$i
  else
    text="IMG_"$i
  fi
  mv $f carpetaNueva/$text".JPG"
  let i=i+1
done
