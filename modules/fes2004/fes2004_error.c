#include <stdio.h>
#include <stdlib.h>

void print_error_1() 
{
  printf("error: wrong run parameters\n use: [executable] [time origine (CNES:0,CTO:1,ESA:2)] [input file] [Nb CPU dispo]\nexit");
  exit(101);

}

void print_error_2() 
{
  printf("a strange error occur, certainly due to the memory allocation in the spectrum initialisation --> exit\n");
  exit(102);

}

void print_error_3(char *message) 
{
  printf("%s\n",message);
  exit(103);

}

void print_error_4(char *message) 
{
  printf("%s\n",message);
  exit(104);

}

void print_error_5(char *message) 
{
  printf("can not open the file : %s , please check this error --> exit\n",message);
  exit(105);

}
