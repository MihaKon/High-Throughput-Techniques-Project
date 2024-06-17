FROM tw_gatk:2024

RUN mkdir /task_data

WORKDIR /task_data

RUN cp /data/NA12878_chr20_1.fastq.gz /task_data/ 
RUN cp /data/NA12878_chr20_2.fastq.gz /task_data/

COPY task_script.sh /task_data/
RUN sed -i 's/\r//' task_script.sh

CMD ["bash", "task_script.sh"]