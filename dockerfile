FROM tw_gatk:2024

RUN mkdir /task_data
WORKDIR /task_data

# CMD [ "sed -i 's/\r//' task_script.sh"]
