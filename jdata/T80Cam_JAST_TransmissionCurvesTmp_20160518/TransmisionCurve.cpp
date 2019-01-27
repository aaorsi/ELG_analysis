#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

//-------------------------------------FUNCIONES ---------------------------------------------------

void generaCurvaTransmision(int valor,char *a[100]);



//-----------------------------------------------------------------------------------------------------

int main()
{
int i;
char *a[100] = {"T80Cam_T80_uJAVA.tab","T80Cam_T80_J0378.tab", "T80Cam_T80_J0395.tab","T80Cam_T80_J0410.tab", "T80Cam_T80_J0430.tab","T80Cam_T80_gSDSS.tab", "T80Cam_T80_J0515.tab", "T80Cam_T80_rSDSS.tab", "T80Cam_T80_J0660.tab", "T80Cam_T80_iSDSS.tab", "T80Cam_T80_J0861.tab", "T80Cam_T80_zSDSS.tab"};
for(i=0;i<12;i++) printf("%s\n",a[i]);
for(i=0;i<12;i++) generaCurvaTransmision(i,a);
}

//FUNCIONES--------------------------------------------------------------------------------------------------------

void generaCurvaTransmision(int valor,char *a[100])
{
FILE *file1;
FILE *fileout;
char caracter;
char nombreFichero[100],nameoutput[100],name[100];
double datoleido,longitudOnda[70000],transisionPlusCCD[70000],central,aux1,aux2,aux3,auxold,oldl,oldt;
double FlongitudOnda[70000],FtransisionPlusCCD[70000];
int i,j,n,k,dim,pos1,pos2;
n=0;
double bin;
bin=25;
char nameaux1[100]="jp";
char nameaux2[100]="_trans.dat";
/*
printf("Give me the file name that I have to do the transmision curve.\n");
scanf("%s",nombreFichero);
*/

file1=fopen(a[valor],"r");

if(file1==NULL) printf("El fichero no existe\n");

else {
fscanf(file1, "%c",&caracter);
if(caracter=='#')
{
    while(caracter!='\n')
    {
        fscanf(file1, "%c",&caracter);
        //printf("%c",caracter);
    }
}

i=0;
fscanf(file1, "%lf",&aux1);
printf("dato leido %lf\n",aux1);
longitudOnda[i]=aux1;
fscanf(file1, "%lf",&aux3);
transisionPlusCCD[i]=aux3;
i++;

fscanf(file1, "%lf",&aux1);
fscanf(file1, "%lf",&aux3);

while(!feof(file1))
{

longitudOnda[i]=aux1;
transisionPlusCCD[i]=aux3;
i++;

fscanf(file1, "%lf",&aux1);
fscanf(file1, "%lf",&aux3);

}
dim=i;
fclose(file1);

}

oldt=transisionPlusCCD[0];

for(i=1;i<dim;i++)
{
if(oldt==0 && transisionPlusCCD[i]>0)
{
pos1=i;
i=dim+10;
}
else oldt=transisionPlusCCD[i];
}

printf("posicion:%d\n",pos1);


oldt=transisionPlusCCD[dim-1];
for(i=dim-2;i>0;i--)
{
    if(oldt==0 && transisionPlusCCD[i]>0)
        {
            pos2=i;
            i=-1;
        }
else oldt=transisionPlusCCD[i];
}

int longitud;
longitud=pos2-pos1;

for(i=0;i<=longitud;i++)
{
FlongitudOnda[i]=longitudOnda[pos1+i];
FtransisionPlusCCD[i]=transisionPlusCCD[pos1+i];

}

//FlongitudOnda[pos1+i]=longitudOnda[pos1+i]+0.10;
//FtransisionPlusCCD[pos1+i]=0.0;

//for(j=0;j<i;j++) printf("longitud de onda=%lf\n",longitudOnda[j]);
sprintf(name,"%s%d",nameaux1,valor+1);
sprintf(nameoutput,"%s%s",name,nameaux2);
fileout=fopen(nameoutput,"w");
fprintf(fileout,"%d\n",longitud+1);
for(j=0;j<=(longitud);j++)
{
fprintf(fileout,"%lf\t%lf\t%lf\t%lf\t%lf\n",FlongitudOnda[j],FtransisionPlusCCD[j],FtransisionPlusCCD[j],FtransisionPlusCCD[j],FtransisionPlusCCD[j]);
}


fclose(fileout);

}
