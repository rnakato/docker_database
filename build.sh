for tag in 2024.04 latest
do
    docker build -t rnakato/database_gpu:$tag --target gpu . #--no-cache
    docker push     rnakato/database_gpu:$tag
    docker build -t rnakato/database:$tag --target normal . #--no-cache
    docker push     rnakato/database:$tag
done
