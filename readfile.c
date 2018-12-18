#include <stdio.h>
#include <stdlib.h>

void main(int argc, char *argv[])
{
  int ch=0;
  int noLines=0;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  float f1,f2;
  char pend1[40];
char pend2[40];

  char *filename= "cities.txt";
  FILE *fp = fopen(filename,"r");

 /*while ( (ch = fgetc(fp)) != EOF )
  {
    if (ch =='\n')
    {
      noLines++;
    //  printf("Line %d\n",noLines);
    }
  } */
  while ((read = getline(&line, &len, fp)) != -1) {
    sscanf(line, "%*s%*c%s%*c%s%*c", pend1,pend2);
    printf("hi  %s  %s\n",pend1,pend2);
      //  printf("%s", line);
  }
  fclose(fp);

  printf("%d\n",noLines);
}
