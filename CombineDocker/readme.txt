docker build . -t combine_splice

docker run -v C:/Splicer/Data:/splice combine_splice /splice ACC

To publish to google:
gcloud auth configure-docker
docker tag combine_splice gcr.io/dockerreg-340916/combine_splice
docker push gcr.io/dockerreg-340916/combine_splice