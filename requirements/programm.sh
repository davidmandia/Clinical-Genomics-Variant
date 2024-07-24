## This Bash script will run the complete database built for both genome version

bash requirements/blast.sh
python requirements/download_SAM_fromS3.py
bash requirements.sh
bash requirements/scripts.sh