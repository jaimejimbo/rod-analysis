for f in *.gif
do
    echo "Converting file $f to mp4..."
    ffmpeg -f gif -i $f $f.mp4
    rm $f
done

