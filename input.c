
int main(int argc,char *argv[]){
    int opt;
    while((opt=getopt(argc,argv,"e:o:t:T:s")) != -1) {
        switch (opt) {
case 'e':
			epsilon=atof(optarg);
			break;
case 'o':
			data=optarg;
			break;
case 't':
            t=atoi(optarg);
            break;
case 'T':
            T=atoi(optarg);
            break;
case 's':
			sflag=1;
			break;
default:
            fprintf(stderr,
"Usage: %s [options]\n"
"Options:\n"
"    -e x    Set noise level to x.\n"
"    -o s    Set output directory to s.\n"
"    -T n    Approximate supremum by taking t less than n.\n"
"    -s      Skip the compile and the generation of data.\n"
"    -t n    Approximate supremum by taking t greater than n.\n",
                argv[0]);
            exit(1);
        }
    }
}
