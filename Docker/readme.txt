docker build . -t splicer

docker run -v C:/Splicer/Validate:/splicer splicer notebooks_GTEX-1117F-0426-SM-5EGHI.Aligned.sortedByCoord.out.patched.md.bam
docker run -v C:/Splicer/:/splicer -w /splicer splicer /splicer/Validate/notebooks_GTEX-1117F-0426-SM-5EGHI.Aligned.sortedByCoord.out.patched.md.bam

To publish to google:
gcloud auth configure-docker
docker tag splicer gcr.io/dockerreg-340916/splicer
docker push gcr.io/dockerreg-340916/splicer