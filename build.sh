for tag in 2022.09 latest #Ensembl106 latest
do
    docker build -t rnakato/database:$tag . #--no-cache
    docker push rnakato/database:$tag
done

for tag in 2022.09 latest #Ensembl106 latest
do
    docker build -f Dockerfile.GPU -t rnakato/database_gpu:$tag . #--no-cache
    docker push rnakato/database_gpu:$tag
done
