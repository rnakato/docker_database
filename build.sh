for tag in 2026.06 latest
do
    docker build -t rnakato/database:$tag --target normal . #--no-cache
#    apptainer build -F /work3/SingularityImages/database.$tag.sif docker-daemon://rnakato/database:$tag
    docker push     rnakato/database:$tag
done

#    docker build -t rnakato/database_gpu:$tag --target gpu . #--no-cache
#    docker push     rnakato/database_gpu:$tag
