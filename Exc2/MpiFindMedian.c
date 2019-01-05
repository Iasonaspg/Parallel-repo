/*

Parallel and Distributed Systems - Assignment 2
Authors: Michail Iason Pavlidis (9015) - Christianou Georgia (8419)
emails: michailpg@ece.auth.gr - gkchristi@ece.auth.gr





The MIT License (MIT)

Copyright (c) 2014

Athanassios Kintsakis
Contact
athanassios.kintsakis@gmail.com
akintsakis@issel.ee.auth.gr

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#define dimension 2

MPI_Status Stat;
// Allagi gia float
void partition (float *array, float* points, int partLength, int elements, float pivot, float **arraysmall, float **arraybig, float** pointsBig, float** pointsSmall, int *endsmall, int *endbig);
float selection(float *array,int number,float *points,int partLength);

/***Kills processes that have no values left in their arrays****/
void removeElement(int *array, int *size, int element)
{
    int i;
    int flag=0;
    for(i=0;i<*size;i++)
    {
        if(flag==1)
            array[i]=array[i+1];
        if(array[i]==element&& flag==0)
        {
            array[i]=array[i+1];
            flag=1;
        }
    }
    *size=*size-1;
}

/***Calculate Lengths and Send them to the corresponding Node***/
void sendLengths(int size,int noProcesses)
{
    int i,partLength;
    if(size%noProcesses!=0)
    {
        int left=size-(size/noProcesses)*noProcesses;  //Split the size in as close to equal as possible parts
        partLength=(size/noProcesses)+1;
        for(i=1;i<left;i++)      //start from 1 because we create the zero one through the main function
            MPI_Send(&partLength,1,MPI_INT,i,1,MPI_COMM_WORLD);
        partLength-=1;
        for(i=left;i<noProcesses;i++)
            MPI_Send(&partLength,1,MPI_INT,i,1,MPI_COMM_WORLD);
    }
    else
    {
        partLength=size/noProcesses;
        for(i=1;i<noProcesses;i++)
            MPI_Send(&partLength,1,MPI_INT,i,1,MPI_COMM_WORLD);
    }
}

/****Swaps two values in an array****/
// Allagi gia float
void swap_values(float *array,float* points,int partLength,int x,int y)
{
    float temp;
    temp=array[x];
    array[x]=array[y];
    array[y]=temp;
    float temp1;
    for (int i=0;i<dimension;i++){
        temp1 = points[i*partLength+x];
        points[i*partLength+x] = points[i*partLength+y];
        points[i*partLength+y] = temp1;
        //printf("Swap: %f with: %f\n", temp1, points[i*partLength+x]);
    }
}

/*****Send random numbers to every node.*****/
// Allagi gia float
void generateNumbers(float *numberPart,int partLength, int cal)
{
    srand((cal+1)*time(NULL));     //Generate number to fill the array, den arkei mono to time() giati trexoun parallhla! Ara to allazoume me to cal pou einai rank!
    int i;
    for(i=0; i<partLength*dimension; i++){
        numberPart[i]= (float)rand();
        if (numberPart[i] > 10000000){
            numberPart[i] = numberPart[i]/348128.24541;
        }
        else if (numberPart[i] > 10000){
            numberPart[i] = numberPart[i]/3481.24541;
        }
        else {
            numberPart[i] = numberPart[i]/34.24541;
        }
        //printf("%f \n",numberPart[i]);
    }
}


/***Validates the stability of the operation****/
// Allagi gia float
void validation(float median,int partLength,int size,float *numberPart,int processId, MPI_Comm sub_comm)
{
    MPI_Bcast(&median,1,MPI_FLOAT,0,sub_comm);
	int countMin=0;
    int countMax=0;
    int countEq=0;
    int sumMax,sumMin,sumEq,i;
    for(i=0;i<partLength;i++)
    {
        if(numberPart[i]>median)
            countMax++;
        else if(numberPart[i]<median)
            countMin++;
        else
            countEq++;
    }
    MPI_Reduce(&countMax,&sumMax,1,MPI_INT,MPI_SUM,0,sub_comm);
    MPI_Reduce(&countMin,&sumMin,1,MPI_INT,MPI_SUM,0,sub_comm);
    MPI_Reduce(&countEq,&sumEq,1,MPI_INT,MPI_SUM,0,sub_comm);
    if(processId==0)
    {
        if((sumMax<=size/2)&&(sumMin<=size/2)){  //Checks if both the lower and higher values occupy less than 50% of the total array.
            printf("VALIDATION PASSED!\n");
        }
        else{
            printf("VALIDATION FAILED!\n");
        }
        printf("Values greater than median: %d\n",sumMax);
        printf("Values equal to median: %d\n",sumEq);
        printf("Values lower than median: %d\n",sumMin);
    }

}

/***Validates the stability of the operation (Single Threaded)****/
// Allagi gia float
void validationST(float median,int size,float *numberPart)
{
	int countMin=0;
    int countMax=0;
    int countEq=0;
    int i;
    for(i=0;i<size;i++)
    {
        if(numberPart[i]>median)
            countMax++;
        else if(numberPart[i]<median)
            countMin++;
        else
            countEq++;
    }
    if((countMax<=size/2)&&(countMin<=size/2)){  //Checks if both the lower and higher values occupy less than 50% of the total array.
        printf("VALIDATION PASSED!\n");
    }
    else{
        printf("VALIDATION FAILED!\n");
    }
	printf("Values greater than median: %d\n",countMax);
    printf("Values equal to median: %d\n",countEq);
    printf("Values lower than median: %d\n",countMin);
}

/****Part executed only by the Master Node****/
// Allagi gia float
float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,float* points, MPI_Comm sub_comm) //MASTER NODE CODE
{
    int elements,i,keepBigSet,sumSets,finalize,randomNode,k;
    float median,pivot,tempPivot;       // Allagi gia float
    int endSmall=0;
    int dropoutFlag=0;
    int endBig=0;
    float *arraySmall,*arrayBig,*arrayToUse;    // Allagi gia float
    float *pointsSmall,*pointsBig,*pointsToUse;
    int* activeNodes;
    int activeSize=noProcesses;
    int stillActive=1;
    int oldSumSets=-1;
    int checkIdentical=0;
    int useNewPivot=0;
    int *pivotArray;
    k=(int)size/2+1; //It is done so in order to find the right median in an even numbered array.
    elements=partLength;
    activeNodes=(int *)malloc(noProcesses*sizeof(int));  //we create the array that contains the active nodes.
    arrayToUse=numberPart;
    pointsToUse = points;
    pivotArray=(int*)malloc(noProcesses*sizeof(int));  //Used for special occasions to gather values different than the pivot.
    for(i=0;i<activeSize;i++)
    {
        activeNodes[i]=i;
    }
    int randomCounter=0;
    int randomCounter2=0;
    struct timeval first, second, lapsed;
    struct timezone tzp;
    gettimeofday(&first, &tzp);
    for(;;)   //Begin the infinite loop until the median is found.
    {
        int counter=0;
        useNewPivot=0;
        if(stillActive==1&&checkIdentical!=0)  //If i still have values in my array and the Sumed Big Set is identical to the previous one, check for identical values.
        {
            for(i=0;i<elements;i++)
            {
                if(pivot==arrayToUse[i])
                    counter++;
                else
                {
                    useNewPivot=1;
                    tempPivot=arrayToUse[i];
                    break;
                }
            }
        }
        if(checkIdentical!=0)
        {
            int useNewPivotMax=0;
	        MPI_Reduce(&useNewPivot,&useNewPivotMax,1,MPI_INT,MPI_MAX,0,sub_comm); //FIRST(OPTIONAL) REDUCE : MAX useNewPivot
            if(useNewPivotMax!=1)    //That means that the only values left are equal to the pivot!
            {
                median=pivot;
                finalize=1;
                MPI_Bcast(&finalize,1,MPI_INT,0,sub_comm); //FIRST(OPTIONAL) BROADCAST : WAIT FOR FINALIZE COMMAND OR NOT
                gettimeofday(&second, &tzp);
                if(first.tv_usec>second.tv_usec)
                {
                    second.tv_usec += 1000000;
                    second.tv_sec--;
                }
                lapsed.tv_usec = second.tv_usec - first.tv_usec;
                lapsed.tv_sec = second.tv_sec - first.tv_sec;
                printf("\n\nTime elapsed: %lu, %lu s\n", lapsed.tv_sec, lapsed.tv_usec);
                validation(median,partLength,size,numberPart,processId,sub_comm);
                //MPI_Finalize();
                free(pivotArray);
                return median;
            }
            else
            {
                finalize=0;
                int useit=0;
                randomCounter2++;
                MPI_Bcast(&finalize,1,MPI_INT,0,sub_comm);
                MPI_Gather(&useNewPivot, 1, MPI_INT, pivotArray, 1, MPI_INT, 0, sub_comm); //Gather every value and chose a node to change the pivot.
                for(i=0;i<activeSize;i++)
                {
                    if(pivotArray[i]==1)
                    {
                        if((randomCounter2>1)&&(randomNode!=activeNodes[i]))  //Check if the same node has already been used in a similar operation.
                        {
                            useit=1;
                            randomNode=activeNodes[i];
                            randomCounter2=0;
                            break;
                        }
                        else if(randomCounter2<2)
                        {
                            useit=1;
                            randomNode=activeNodes[i];
                            break;
                        }
                    }
                }
                if(useit!=0)
                    useNewPivot=1;
                else
                    useNewPivot=0;
            }
        }
        if(useNewPivot!=0){
            MPI_Bcast(&randomNode,1,MPI_INT,0,sub_comm);  //THIRD(OPTIONAL) BROADCAST : BROADCAST THE SPECIAL NODE
        }
        if(useNewPivot==0)  //if we didnt choose a special Node, choose the node that will pick the pivot in a clockwise manner. Only selects one of the active nodes.
        {
            if(randomCounter>=activeSize)
                randomCounter=0; //Fail safe
            randomNode=activeNodes[randomCounter];
            randomCounter++;			//Increase the counter
            MPI_Bcast(&randomNode,1,MPI_INT,0,sub_comm);   //FIRST BROADCAST : SENDING randomnode, who will chose
        }
        if(randomNode==processId)  //If i am to choose the pivot.....
	    {
            if(useNewPivot==0)
            {
                srand(time(NULL));
                pivot=arrayToUse[rand() % elements];
                MPI_Bcast(&pivot,1,MPI_FLOAT,0,sub_comm); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao           // Allagi gia float
	        }
            else
            {
                MPI_Bcast(&tempPivot,1,MPI_FLOAT,0,sub_comm); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao       // Allagi gia float
                pivot=tempPivot;
            }
        }
        else //If not.. wait for the pivot to be received.
            MPI_Bcast(&pivot,1,MPI_FLOAT,randomNode,sub_comm);  // SECOND BROADCAST : RECEIVING PIVOT         // Allagi gia float
        if(stillActive==1)  //If i still have values in my array.. proceed
        {
            partition(arrayToUse,pointsToUse,partLength,elements,pivot,&arraySmall,&arrayBig,&pointsBig,&pointsSmall,&endSmall,&endBig);  //I partition my array  
            // endsmall=number of elements in small array, it may be 0
            // endbig=number of elements in big array, it may be 0
            //arraysmall = Points to the position of the small array.NULL if the array is empty
            //Same for arraybig
        }
        else  //If i'm not active endBig/endSmall has zero value.
        {
            endBig=0;
            endSmall=0;
        }
        sumSets=0;
	    //We add the bigSet Values to decide if we keep the small or the big array
	    MPI_Reduce(&endBig,&sumSets,1,MPI_INT,MPI_SUM,0,sub_comm);  //FIRST REDUCE : SUM OF BIG
        MPI_Bcast(&sumSets,1,MPI_INT,0,sub_comm);
        if(oldSumSets==sumSets)
            checkIdentical=1;
        else
        {
            oldSumSets=sumSets;
            checkIdentical=0;
        }
	    //hmetabliti keepBigSet 0 h 1 einai boolean k me autin enimerwnw ton lao ti na kratisei to bigset h to smallset
	    if(sumSets>k)   //an to sumofbigsets > k tote krataw to big SET
	    {
            keepBigSet=1; //to dilwnw auto gt meta tha to steilw se olous
            if(endBig==0)
                dropoutFlag=1; //wraia.. edw an dw oti to bigset mou einai 0.. alla prepei na kratisw to bigset sikwnw auti ti simaia pou simainei tha ginw inactive ligo pio katw tha to deis
            else
            {
                arrayToUse=arrayBig; //thetw ton neo pinaka na einai o big
                pointsToUse = pointsBig;
                elements=endBig; //thetw arithmo stoixeiwn iso me tou big
            }
	    }
	    else if(sumSets<k) //antistoixa an to sumofbigsets < k tote krataw to small set
	    {
		    keepBigSet=0;
		    k=k-sumSets;
		    if(endSmall==0)
                dropoutFlag=1; //antistoixa koitaw an tha ginw inactive..
		    else
		    {
		    	arrayToUse=arraySmall; //dinw times..
                pointsToUse = pointsSmall;
		    	elements=endSmall;
		    }
	    }
	    else  //edw simainei k=sumofbigsetes ara briskw pivot k telos
	    {
		    median=pivot;
		    finalize=1; //dilwnw finalaize =1
            MPI_Bcast(&finalize,1,MPI_INT,0,sub_comm); //to stelnw se olous, oi opoioi an laboun finalize =1 tote kaloun MPI finalize k telos
		    gettimeofday(&second, &tzp);
            if(first.tv_usec>second.tv_usec)
            {
                second.tv_usec += 1000000;
                second.tv_sec--;
            }
            lapsed.tv_usec = second.tv_usec - first.tv_usec;
            lapsed.tv_sec = second.tv_sec - first.tv_sec;
            printf("\n\nTime elapsed: %lu, %lu s\n", lapsed.tv_sec, lapsed.tv_usec);
		    validation(median,partLength,size,numberPart,processId,sub_comm);
            //MPI_Finalize();
            free(pivotArray);
            return median;
        }
        finalize=0; //an den exw mpei sta if den exw steilei timi gia finalize.. oi alloi omws perimenoun na laboun kati, stelnw loipon to 0 pou simainei sunexizoume
        MPI_Bcast(&finalize,1,MPI_INT,0,sub_comm);	//SECOND BROADCAST : WAIT FOR FINALIZE COMMAND OR NOT
        //edw tous stelnw to keepbigset gia na doun ti tha dialeksoun
        MPI_Bcast(&keepBigSet,1,MPI_INT,0,sub_comm);    //THIRD BROADCAST: SEND keepBigset boolean
        if(dropoutFlag==1 && stillActive==1) //edw sumfwna me to dropoutflag pou orisame prin an einai 1 kalw tin sinartisi pou me petaei apo ton pinaka. episis koitaw na eimai active gt an me exei idi petaksei se proigoumeni epanalispi tote den xreiazetai na me ksanapetaksei
        {
            stillActive=0;
            removeElement(activeNodes, &activeSize, 0);
        }
        int flag;
        //edw perimenw na akousw apo ton kathena an sunexizei active h oxi.. an oxi ton petaw.. an einai idi inactive apo prin stelnei kati allo (oxi 1)k den ton ksanapetaw
        for(i=0;i<activeSize;i++)
        {
            if(activeNodes[i]!=0)
            {
                MPI_Recv(&flag,1,MPI_INT,activeNodes[i],1,sub_comm,&Stat);  //FIRST RECEIVE : RECEIVE active or not
                if(flag==1)
                    removeElement(activeNodes, &activeSize, activeNodes[i]);
            }
        }
    }
}


/***Executed only by Slave nodes!!*****/
// Allagi gia float
void slavePart(int processId,int partLength,float *numberPart,int size,float* points, MPI_Comm sub_comm)  //code here is for the cheap slaves :P
{
	int dropoutflag,elements,i,sumSets,finalize,keepBigSet,randomNode;
    float pivot,tempPivot;          // Allagi gia float
    int endSmall=0;
    int endBig=0;
    float *arraySmall,*arrayBig,*arrayToUse;        // Allagi gia float
    float *pointsSmall,*pointsBig,*pointsToUse;
	arrayToUse=numberPart;
    pointsToUse = points;
	elements=partLength;
	int stillActive=1;
	int *pivotArray;
    int oldSumSets=-1;
    int checkIdentical=0;
    int useNewPivot;
    for(;;)
	{
        finalize=0;
        int counter=0;
        useNewPivot=0;
        if(stillActive==1&&checkIdentical!=0)  //If i still have values in my array..   If the Sumed Big Set is identical to the previous one, check for identical values.
        {
            for(i=0;i<elements;i++)
            {
                if(pivot==arrayToUse[i])
                    counter++;
                else
                {
                    useNewPivot=1;
                    tempPivot=arrayToUse[i];
                    break;
                }
            }
        }
        if(checkIdentical!=0)
        {
            int useNewPivotMax=0;
            MPI_Reduce(&useNewPivot,&useNewPivotMax,1,MPI_INT,MPI_MAX,0,sub_comm);
            MPI_Bcast(&finalize,1,MPI_INT,0,sub_comm);//an o master apo to keepbigset k apo to count apofasisei oti teleiwsame mou stelnei 1, alliws 0 sunexizoume
            if(finalize==1)
            {
                float median=0.0;
                validation(median,partLength,size,numberPart,processId,sub_comm);
                //MPI_Finalize();
                return ;
            }
            else
            {
                MPI_Gather(&useNewPivot, 1, MPI_INT, pivotArray, 1, MPI_INT, 0, sub_comm);
            }
        }
        MPI_Bcast(&randomNode,1,MPI_INT,0,sub_comm); //FIRST BROAD CAST : RECEIVING RANDOM NODE, perimenw na dw poios einaito done
        if(randomNode!=processId) //means I am not the one to chose pivot.. so I wait to receive the pivot
            MPI_Bcast(&pivot,1,MPI_FLOAT,randomNode,sub_comm);	//SECOND BROADCAST : RECEIVING PIVOT
        else if(randomNode==processId) //I am choosing suckers
        {
            if(useNewPivot==0)
            {
                srand(time(NULL));
                pivot=arrayToUse[rand() % elements];
                MPI_Bcast(&pivot,1,MPI_FLOAT,processId,sub_comm); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao       // Allagi gia float
            }
            else
            {
                MPI_Bcast(&tempPivot,1,MPI_FLOAT,processId,sub_comm); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao       // Allagi gia float
                pivot=tempPivot;
            }
        }
        if(stillActive==1)   //an eksakolouthw na eimai active, trexw tin partition.. k to count kommati to opio eimape kapou exei problima
        {
            partition(arrayToUse,pointsToUse,partLength,elements,pivot,&arraySmall,&arrayBig,&pointsBig,&pointsSmall,&endSmall,&endBig);
        }
        else
        {
            endBig=0;
            endSmall=0;
        }
        //an eimai inactive stelnw endbig=0 gia to bigset pou den epireazei
        sumSets=0;
        MPI_Reduce(&endBig,&sumSets,1,MPI_INT,MPI_SUM,0,sub_comm); //FIRST REDUCE : SUM OF BIG, stelnw ola ta bigset gia na athroistoun sotn master
        MPI_Bcast(&sumSets,1,MPI_INT,0,sub_comm);
        if(oldSumSets==sumSets)
            checkIdentical=1;
        else
        {
            oldSumSets=sumSets;
            checkIdentical=0;
        }
        MPI_Bcast(&finalize,1,MPI_INT,0,sub_comm);//an o master apo to keepbigset k apo to count apofasisei oti teleiwsame mou stelnei 1, alliws 0 sunexizoume
        if(finalize==1)
        {
            float median=0.0;
            validation(median,partLength,size,numberPart,processId,sub_comm);
            //MPI_Finalize();
            return ;
        }
        MPI_Bcast(&keepBigSet,1,MPI_INT,0,sub_comm);//THIRD BROADCAST: Receive keepBigset boolean, edw lambanw an krataw to mikro i megalo set.
            //afou elaba ton keepbigset an eimai active krataw enan apo tous duo pinake small h big.. alliws den kanw tpt
            //edw antistoixa allazw tous pointers, k eksetazw an exw meinei xwris stoixeia tin opoia periptwsi sikwnw to dropoutflag k pio katw tha dilwsw na ginw inactive
        if(stillActive==1)
        {
            if(keepBigSet==1)
            {
                if(endBig==0)
                    dropoutflag=1;
                else
                {
                    arrayToUse=arrayBig;
                    pointsToUse = pointsBig;
                    elements=endBig;
                }
            }
            else if(keepBigSet==0)
            {
                if(endSmall==0)
                    dropoutflag=1;
                else
                {
                    arrayToUse=arraySmall;
                    pointsToUse = pointsSmall;
                    elements=endSmall;
                }
            }
        }
        //edw einai ligo periploka grammeno, isws exei perita mesa alla, an eimai active k thelw na ginw inactive einai i prwti periptwsi, h deuteri einai eimai inactive hdh k i triti einai sunexizw dunamika
        if(dropoutflag==1 && stillActive==1)
        {
            MPI_Send(&dropoutflag,1,MPI_INT,0,1,sub_comm); //FIRST SEND : send active or not;
            stillActive=0;
        }
        else if(stillActive==0)
        {
            dropoutflag=-1;
            MPI_Send(&dropoutflag,1,MPI_INT,0,1,sub_comm); //FIRST SEND : send active or not;
        }
        else
        {
            dropoutflag=0;
            MPI_Send(&dropoutflag,1,MPI_INT,0,1,sub_comm); //FIRST SEND : send active or not;
        }
    }
}


/* My functions */

// Selects vantage point and bcast it to your slaves
void selectVP(int processId, float* points, int partLength, int k, float* vp){
    srand((processId+k+1)*time(0));
    int randPo = rand()%partLength;
    // int randPo = 1;
    // printf("Random Position of Vantage Point: %d\n", randPo);
    for (int i=0; i<dimension;i++){
        vp[i] = points[i*partLength+randPo];
    }
}

// Calculates distance from vantage point supporting dynamic dimension size
void calcDistance(float* points, int partLength, float* vp, float* distances){
    for (int i=0; i<partLength; i++){
        distances[i] = 0;
    }
    
    for (int i=0; i<partLength; i++){
        float sum = 0;
        for (int j=0; j<dimension; j++){
            sum += powf(points[j*partLength+i]-vp[j],2);
        }
        distances[i] = powf(sum,0.5);
        //printf("Inside function: %f\n",distances[i]);
    }
}

// Checks for lower, greater and equal values to median, master gathers them. It also rearranges the array in order to be [low:equ:gre]
void checkSums(float* distances,float* points, int partLength, float median, int* gtMedian, int* eqMedian, int* ltMedian, MPI_Comm sub_comm){
    int countMin=0;
    int countMax=0;
    int countEq=0;
    int *posEq = (int*)malloc(partLength*sizeof(int));
    for(int i=0;i<partLength;i++)
    {
        if(distances[i]>median)
            countMax++;
        else if(distances[i]<median)
            countMin++;
        else{
            // posEq = (int*)realloc(posEq,countEq);
            posEq[countEq] = i;
            countEq++;
            if (countEq > 1){
                //printf("Bad dataset values, multiple values equal to median. Program is going to terminate...\n");
                // exit(1);
            }
        }
    }
    for (int i=0;i<countEq;i++){
        swap_values(distances,points,partLength,posEq[i],countMin+i);
        // printf("PosEq: %d\n",posEq[i]);
    }
    free(posEq);
    countMin += countEq;
    MPI_Gather(&countMin,1,MPI_INT,ltMedian,1,MPI_INT,0,sub_comm);  //Gather those values to master
    MPI_Gather(&countEq,1,MPI_INT,eqMedian,1,MPI_INT,0,sub_comm);
    MPI_Gather(&countMax,1,MPI_INT,gtMedian,1,MPI_INT,0,sub_comm);
}

// Serial implmentation for checkSum
void serialCheckSums(float* distances, float* points, int partLength, float median, int* countMin, int* countMax){
    int countEq=0;
    int *posEq = (int*)malloc(partLength*sizeof(int));
    for(int i=0;i<partLength;i++)
    {
        if(distances[i]>median)
            (*countMax)++;
        else if(distances[i]<median)
            (*countMin)++;
        else{
            // posEq = (int*)realloc(posEq,countEq);
            posEq[countEq] = i;
            countEq++;
            if (countEq > 1){
                //printf("Bad dataset values, multiple values equal to median. Program is going to terminate...\n");
                // exit(1);
            }
        }
    }
    for (int i=0;i<countEq;i++){
        swap_values(distances,points,partLength,posEq[i],(*countMin)+i);
        // printf("PosEq: %d\n",posEq[i]);
    }
    if((*countMax<=partLength/2)&&(*countMin<=partLength/2)){  //Checks if both the lower and higher values occupy less than 50% of the total array.
        // printf("Single thread validation PASSED!\n");
    }
    else{
        // printf("Single thread validation FAILED!\n");
    }
    *countMin += countEq;
    free(posEq);
}


// Validates the distribution of the points
void distrValidation(int sub_id, int sub_size, float* distances, float median, int partLength, MPI_Comm sub_comm){
    int flag = 1;
    // int sumFlags = 0;
    
    if (sub_id < sub_size/2){
        for (int i=0; i<partLength; i++){
            if (distances[i] > median){
                printf("Validation of process %d after distribute has failed\n",sub_id);
                flag = 0;
                break;
            }
        }
    }
    else{
        for (int i=0; i<partLength; i++){
            if (distances[i] <= median){
                printf("Validation of process %d after distribute has failed\n",sub_id);
                flag = 0;
                break;
            }
        }
    }

    if (flag){
        printf("Validation of process %d after distribute has passed\n",sub_id);
    }

    // MPI_Reduce(&flag,&sumFlags,1,MPI_INT,MPI_SUM,0,sub_comm);
    // if(sub_id == 0)
    // {
    //     if(sumFlags == sub_size){
    //         printf("All the processes of this communicator have passed the validation!\n");
    //     }
    //     else{
    //         printf("Not all the processes of this communicator have passed the validation!\n");
    //     }  
    // }
}

void Distribute(float* points, int* ltMedian, int* gtMedian, int sub_id, int sub_size, int partLength, MPI_Comm sub_comm){

    int destination;
    // O master theloume na kathorizei tis epikoinonies
    if (sub_id == 0){
        int masterLength = gtMedian[0];     //Thelei na dosei kai na parei tosa osa einai ta megala tou
        int sendStart = ltMedian[0];        //Ta megala ksekinane ston pinaka ekei pou telionoun ta mikrotera/isa
        int recvStart = 0;                  //Apothikeuo apo tin arxi tou tempΒuf
        float *tempBuf = (float*)malloc(masterLength*dimension*sizeof(float));    //Pinakas gia ta recv tou patera
        for (int i=0;i<(sub_size/2);i++){
            int sendCount = gtMedian[i];    //Oi protoi misoi theloun na steiloun ta megala
                        
            for (int j=(sub_size/2);j<sub_size;j++){
                
                if (sendCount == 0){
                    break;              //An o i den thelei na dosei, des gia ton epomeno
                }

                if (ltMedian[j]>0){
                    int toSend;     //Posa mporei na steilei ston j
                    if (ltMedian[j] >= sendCount){
                        toSend = sendCount;
                    }
                    else{
                        toSend = ltMedian[j];
                    }

                    MPI_Send(&i,1,MPI_INT,j,j,sub_comm);    //Leme ston j oti tha kanei antallagi me ton i, toSend stoixeia
                    MPI_Send(&toSend,1,MPI_INT,j,j,sub_comm);

                    if (i != 0){
                        MPI_Send(&j,1,MPI_INT,i,i,sub_comm);    //Leme ston i oti tha kanei antallagi me ton j, toSend plithos
                        MPI_Send(&toSend,1,MPI_INT,i,i,sub_comm);
                    }
                    else{
                        destination = j;    //Ksekiname epi topou tin antallagi tou master
                        
                        for (int l=0;l<dimension;l++){
                            MPI_Sendrecv(&points[partLength*l+sendStart], toSend, MPI_FLOAT, destination, sub_id,&tempBuf[masterLength*l + recvStart],toSend,MPI_FLOAT,destination,destination,sub_comm,MPI_STATUS_IGNORE);
                        }

                        // for (int i=0; i<toSend; i++){
                        //     for (int l=0; l<dimension; l++){
                        //         printf("Master receive: %.7f\n",(&tempBuf[masterLength*l + recvStart])[i]);
                        //     }
                        // }
    
                        // for (int i=0; i<toSend; i++){
                        //     for (int l=0; l<dimension; l++){
                        //         printf("Master sending: %.7f\n",(&points[partLength*l + sendStart])[i]);
                        //     }
                        // }

                        // printf("Master exchanged %d points with process %d\n",toSend,destination);

                    }

                    sendStart += toSend;    //Steilame ta prota toSend stoixeia
                    recvStart += toSend;    //Lavame ta prota toSend stoixeia
                    sendCount -= toSend;    //Pleon thelei na steilei ligotera
                    ltMedian[j] -= toSend;    //Pleon o j thelei na parei ligotera  
                }
            }
        }

        for (int i=0; i<masterLength; i++){
            for (int j=0; j<dimension; j++){
                (&points[j*partLength + ltMedian[0]])[i] = tempBuf[j*masterLength + i];     //Ensomatonoume ta stoixeia pou lavame
            }
        }


    }
    else{

        if(sub_id < sub_size/2){ 
            
            int sendCount = gtMedian[sub_id];    //Oi prwtoi stelnoun ta megala
            int slaveLength = gtMedian[sub_id];     //Posa thelei na parei kai na dosei sunolika
            int sendStart = ltMedian[sub_id];       //Ta megala ksekinane ston pinaka ekei pou telionoun ta mikrotera/isa
            int recvStart = 0;
            int toSend;
            float *tempBuf = (float*)malloc(slaveLength*dimension*sizeof(float));

            while (sendCount){
                
                MPI_Recv(&destination,1,MPI_INT,0,sub_id,sub_comm,MPI_STATUS_IGNORE);   //O master leei me poion tha ginei antallagi kai to plithos
                MPI_Recv(&toSend,1,MPI_INT,0,sub_id,sub_comm,MPI_STATUS_IGNORE);
                
                for (int l=0;l<dimension;l++){
                    MPI_Sendrecv(&points[partLength*l+sendStart], toSend, MPI_FLOAT, destination, sub_id,&tempBuf[slaveLength*l + recvStart],toSend,MPI_FLOAT,destination,destination,sub_comm,MPI_STATUS_IGNORE);
                }

                // for (int i=0; i<toSend; i++){
                //     for (int l=0; l<dimension; l++){
                //         printf("First half slave receive: %.7f\n",(&tempBuf[slaveLength*l + recvStart])[i]);
                //     }
                // }
    
                // for (int i=0; i<toSend; i++){
                //     for (int l=0; l<dimension; l++){
                //         printf("First half slave sending: %.7f\n",(&points[partLength*l + sendStart])[i]);
                //     }
                // }

                sendStart += toSend;     //Steilame ta prota toSend stoixeia
                recvStart += toSend;     //Lavame ta prota toSend stoixeia
                sendCount -= toSend;     //Afairoume auta pou steilame

                // printf("Process %d exchanged %d points with process %d. Remaining points: %d\n",sub_id,toSend,destination,sendCount);
            }

            for (int i=0; i<slaveLength; i++){
                for (int j=0; j<dimension; j++){
                    (&points[j*partLength + ltMedian[sub_id]])[i] = tempBuf[j*slaveLength + i];     //Ensomatonoume ta stoixeia pou lavame
                }
            }
            
        }
        else{
            
            int sendCount = ltMedian[sub_id];      //Oi upoloipoi stelnoun ta mikra
            int slaveLength = ltMedian[sub_id];
            int start = 0;                          //Ta mikrotera/isa ksekinan apo tin arxi tou pinaka
            int toSend;
            float *tempBuf = (float*)malloc(slaveLength*dimension*sizeof(float));

            while (sendCount){

                                
                MPI_Recv(&destination,1,MPI_INT,0,sub_id,sub_comm,MPI_STATUS_IGNORE);
                MPI_Recv(&toSend,1,MPI_INT,0,sub_id,sub_comm,MPI_STATUS_IGNORE);
                
                for (int l=0;l<dimension;l++){
                    MPI_Sendrecv(&points[partLength*l+start], toSend, MPI_FLOAT, destination, sub_id,&tempBuf[slaveLength*l + start],toSend,MPI_FLOAT,destination,destination,sub_comm,MPI_STATUS_IGNORE);
                }

                // for (int i=0; i<toSend; i++){
                //     for (int l=0; l<dimension; l++){
                //         printf("Slave sending: %.7f\n",(&points[partLength*l + start])[i]);
                //     }
                // }

                // for (int i=0; i<toSend; i++){
                //     for (int l=0; l<dimension; l++){
                //         printf("Slave receive: %.7f\n",(&tempBuf[slaveLength*l + start])[i]);
                //     }
                // }
                             
                start += toSend;        //Steilame/lavame ta prota toSend stoixeia
                sendCount -= toSend;    //Afairoume osa steilame

                // printf("Process %d exchanged %d points with process %d. Remaining points: %d\n",sub_id,toSend,destination,sendCount);
            }

            for (int i=0; i<slaveLength; i++){
                for (int j=0; j<dimension; j++){
                    points[j*partLength+i] = tempBuf[j*slaveLength + i];    //Ensomatonoume ta stoixeia pou lavame
                }
            }
        }
    }
}

// Swap Points based on "sorted" distances after partition, Too slow for big size arrays
void swapPoints(float* points, float* distances, float* distances_bk, int partLength){
    int* aux = (int*)malloc(partLength*sizeof(int));
    int* flags = (int*)malloc(partLength*sizeof(int));
    for (int i=0; i<partLength; i++){
        flags[i] = 0;
    }
    for (int i=0; i<partLength; i++){
        for (int j=0; j<partLength; j++){
            if ((distances[j] == distances_bk[i]) && (flags[j] == 0)){
                aux[i] = j;
                flags[j] = 1;
                break;
            }
        }
    }
    float* pointsAux = (float*)malloc(partLength*dimension*sizeof(float));
    for (int i=0; i<partLength; i++){
        for (int j=0; j<dimension; j++){
            pointsAux[j*partLength + aux[i]] = points[j*partLength + i];
        }
    }
    memcpy(points,pointsAux,partLength*dimension*sizeof(float));
    free(aux);
    free(flags);
    free(pointsAux);
}

// Creates local VP trees
void serialVP(float *points, float *localVPs, float *localMedians, int size, int subSize, int index, int processId) 
{
    if (subSize > 2) 
    {

        if (index>size){
            // printf("Index: %d\n",index);
            return;
        } 
        int countMin = 0;
        int countMax = 0;
        float* vp = (float*)malloc(dimension*sizeof(float));
        selectVP(processId,points,subSize,index,vp);
        // printf("vantage point inside function: (%f,%f)\n",vp[0],vp[1]);

        for (int i=0; i<dimension; i++){
            localVPs[i*size + index] = vp[i];   //store vp in the "tree"
        }
        
        float* distances = (float*)malloc(subSize*sizeof(float));      
        calcDistance(points,subSize,vp,distances);

        float median = selection(distances,subSize,points,subSize);
        serialCheckSums(distances,points,subSize,median,&countMin,&countMax);
        
        // printf("Median: %f\n",median);

        for (int i=0; i<subSize; i++){
            // printf("distance: %f\n",distances[i]);
        }
        // for (int i=0; i<size; i++){
        //     printf("points after median: %f\n",points[i]);
        // }
        localMedians[index]=median;     //store median in the "tree"

            
        free(distances);
        free(vp);

        serialVP(points,localVPs,localMedians,size,countMin,index*2,processId);     //left child
        
        serialVP(&points[countMin],localVPs,localMedians,size,countMax,index*2+1,processId);    //right child
    }
}


/*****MAIN!!!!!!!!!!*****/
// Allagi gia float
int main (int argc, char **argv)
{
    int processId,noProcesses,size,partLength;
    float median;
    float *points;
    float* localVPs;      //VPs for local trees
    float* localMedians;  //Medians for local trees
    size= (1<<atoi(argv[1]));
    MPI_Init (&argc, &argv);	/* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &processId);	/* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &noProcesses);	/* get number of processes */
    float *dataX = NULL;
    float *dataY = NULL;
    double start_time,end_time;
    


    if(processId==0)
    {
        printf("size: %d processes: %d\n",size,noProcesses);
        
        if(noProcesses>1)
        {
            if(size%noProcesses==0)
                partLength=(size/noProcesses);
            else
                partLength=(size/noProcesses)+1;
            
            sendLengths(size,noProcesses);
            points=(float*)malloc(dimension*partLength*sizeof(float));    //Desmeuoyme xwro gia ton pinaka otan einai parallhlo
            // generateNumbers(points,partLength,processId);


            // Read data set
            FILE * fp1;
            FILE * fp2;
            char *line1 = NULL;
            char *line2 = NULL;
            size_t len1 = 0;
            size_t len2 = 0;
            dataX = (float*)malloc((size)*sizeof(float));
            dataY = (float*)malloc((size)*sizeof(float));
        
            fp1 = fopen("./Dim1.txt", "r");
            fp2 = fopen("./Dim2.txt", "r");
            if ((fp1 == NULL) || (fp2 == NULL)){
                exit(EXIT_FAILURE);
            }
                
            for (int k=0; k<size; k++) {
                getline(&line1, &len1, fp1);
                dataX[k] = atof(line1);
                getline(&line2, &len2, fp2);
                dataY[k] = atof(line2);
                // printf("%s\n", line);
                //printf("FloatX: %f\n",dataX[k]);
                //printf("FloatY: %f\n",dataY[k]);
            }
        
            fclose(fp1);
            fclose(fp2);
            
            if (line1){
                free(line1);
            }
            if (line2){
                free(line2);
            }
            // End of reading data set      

            // Scatter the global array to slaves
            MPI_Scatter(dataX, partLength, MPI_FLOAT, points, partLength, MPI_FLOAT, 0, MPI_COMM_WORLD);
            MPI_Scatter(dataY, partLength, MPI_FLOAT, &points[partLength], partLength, MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
        else
        {
            points = (float*)malloc(dimension*size*sizeof(float));  //Desmeuoyme xwro gia ton pinaka otan einai seiriako
            localVPs = (float*)malloc(dimension*size*sizeof(float));    //Desmeuoyme xwro gia ta Vantage points kai tous medians.
            localMedians = (float*)malloc(size*sizeof(float));      //Zitame mia thesi parapanw apo auto pou xreiazomaste giati paraleipoume tin thesi 0
            // generateNumbers(points,size,processId); //Arxikopoioume ton pinaka me times
            
            // Read data set
            FILE * fp1;
            FILE * fp2;
            char *line1 = NULL;
            char *line2 = NULL;
            size_t len1 = 0;
            size_t len2 = 0;
                   
            fp1 = fopen("./Dim1.txt", "r");
            fp2 = fopen("./Dim2.txt", "r");
            if ((fp1 == NULL) || (fp2 == NULL)){
                exit(EXIT_FAILURE);
            }
                
            for (int k=0; k<size; k++) {
                getline(&line1, &len1, fp1);
                points[k] = atof(line1);
                getline(&line2, &len2, fp2);
                points[size + k] = atof(line2);
                // printf("%s\n", line);
                //printf("FloatX: %f\n",dataX[k]);
                //printf("FloatY: %f\n",dataY[k]);
            }
        
            fclose(fp1);
            fclose(fp2);
            
            if (line1){
                free(line1);
            }
            if (line2){
                free(line2);
            }
            // End of reading data set

            for (int i=0; i<dimension*size; i++){
                //printf("Points from serial: %f\n",points[i]);
            }
            printf("Sequential tree creation has begun\n");
            start_time = MPI_Wtime();
            serialVP(points,localVPs,localMedians,size,size,1,1);
            end_time = MPI_Wtime();
            printf("Sequential tree creation has been completed\n");

            double Worktime = end_time - start_time;

            if (size > 20){
                printf("First 10 Vantage Points of the local tree:\n");
                for (int i=1; i<11; i++){
                     printf("Local Vantage Point from local tree: (%f,%f)\n",localVPs[i],localVPs[size + i]);
                }
                printf("First 10 Medians of the local tree:\n");
                for (int i=1; i<11; i++){
                    printf("Local Median from local tree: %f\n",localMedians[i]);
                }
            }
            else{
                for (int i=1; i<size/2; i++){
                    printf("Local Vantage Point from local tree: (%f,%f)\n",localVPs[i],localVPs[size + i]);
                }
                for (int i=1; i<size/2; i++){
                    printf("Local Medians from serial: %f\n",localMedians[i]);
                }
            }
                       
            printf("Sequential Time elapsed: %f, s\n", Worktime);
            free(points);
            free(localVPs);
            free(localMedians);
            MPI_Finalize();
            return 0;
        }
    }
    else
    {
        MPI_Recv(&partLength,1,MPI_INT,0,1,MPI_COMM_WORLD,&Stat);   //kodikas pou trexei se kathe paidi
        points=(float*)malloc(dimension*partLength*sizeof(float));        //Oi metavlhtes einai topikes se kathe proc
        MPI_Scatter(dataX, partLength, MPI_FLOAT, points, partLength, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Scatter(dataY, partLength, MPI_FLOAT, &points[partLength], partLength, MPI_FLOAT, 0, MPI_COMM_WORLD);
        // generateNumbers(points,partLength,processId);
    }

    //End of data scattering
        
    int sub_id, sub_size,masters_id, masters_size;
    int threshold = log2(noProcesses); //depth of the tree
    float* sharedVPs = (float*)malloc((noProcesses)*dimension*sizeof(float));   //Pinakes gia ta vantage points kai tous medians tou koinou dentrou
    float* sharedMedians = (float*)malloc((noProcesses)*sizeof(float));     //Zitame mia thesi parapanw giati paraleipoume ti thesi 0
    localVPs = (float*)malloc(partLength*dimension*sizeof(float));          ////Pinakes gia ta vantage points kai tous medians tou local dentrou
    localMedians = (float*)malloc(partLength*sizeof(float));
    int start = 0;  //counter gia tin apothikeusi sto dentro
 
    MPI_Barrier(MPI_COMM_WORLD);
    if (processId == 0){
        printf("Tree construction has begun!\n");
        start_time = MPI_Wtime();
    }
    for (int k=0; k<threshold; k++){
        
        int color = processId/(noProcesses/(1<<k)); // Determine color based on level
        
        // Split the communicator based on the color and use the original rank for ordering
        MPI_Comm sub_comm,masters_comm;
        MPI_Comm_split(MPI_COMM_WORLD, color, processId, &sub_comm);
        
        MPI_Comm_rank(sub_comm, &sub_id);
        MPI_Comm_size(sub_comm, &sub_size);
        int *gtMedian = (int*)malloc(sub_size*sizeof(int)); 
        int *eqMedian = (int*)malloc(sub_size*sizeof(int));     //Pinakes pou tha deixnoun posa low/equ/gre tou median exei to kathe processe
        int *ltMedian = (int*)malloc(sub_size*sizeof(int));
        
        if (sub_id == 0){
            MPI_Comm_split(MPI_COMM_WORLD, 0, processId, &masters_comm);        // Ftiaxnoume communicator mono gia tous masters

            // MPI_Comm_rank(masters_comm, &masters_id);

            // MPI_Comm_size(masters_comm, &masters_size);
        }
        else{
            MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, processId, &masters_comm);    // Do not include masters in the communicator
        }

       
        if(sub_id == 0)
        {
            
            // printf("MPI_COMM_WORLD RANK/SIZE: %d/%d \t MPI_CAPTAINS_COMM RANK/SIZE: %d/%d\n", processId, noProcesses, masters_id, masters_size);

            float* vp = (float*)malloc(dimension*sizeof(float));
            selectVP(processId,points,partLength,k,vp);
            MPI_Bcast(vp,dimension,MPI_FLOAT,0,sub_comm);
            for (int j=0; j<dimension; j++){
                MPI_Gather(&vp[j], 1, MPI_FLOAT, &sharedVPs[j*(noProcesses-1) + start + 1], 1, MPI_FLOAT, 0, masters_comm);  // PID=0 gathers the VPs from masters
            }
            printf("Vantage point at level-%d: (%.7f,%.7f)\n",k, vp[0], vp[1]);
            float* distances = (float*)malloc(partLength*sizeof(float));
            //float* distances_bk = (float*)malloc(partLength*sizeof(float));
            calcDistance(points, partLength, vp, distances);
            //memcpy(distances_bk,distances,partLength*sizeof(float));
            // for (int p=0; p<partLength; p++){
            //     printf("Point from master before median: (%.7f,%.7f)\n",points[p],points[partLength + p]);
            // }
            // for (int p=0; p<partLength; p++){
            //     printf("Distance from master before median: %.7f\n",distances[p]);
            // }
            median=masterPart(sub_size,sub_id,size/(1<<k),partLength,distances,points,sub_comm);        // Vriskoume median
            MPI_Bcast(&median,1,MPI_FLOAT,0,sub_comm);                                                  // Bcast him to every slave in the team
            MPI_Gather(&median, 1, MPI_FLOAT, &sharedMedians[start + 1], 1, MPI_FLOAT, 0, masters_comm);  // PID=0 gathers the medians from masters
            start += (1<<k);                                                                               // Se kathe epipedo exoume 2^k medians/VPs
            printf("Median: %.7f\n\n",median);  // Allagi gia float
            checkSums(distances,points,partLength,median,gtMedian,eqMedian,ltMedian,sub_comm);          // Elegxoume ek neou ta sums kai kanoume rearrange an xreiazetai
            //swapPoints(points,distances,distances_bk,partLength,dimension);
            // for (int p=0; p<partLength; p++){
            //     printf("Point from master before distribute: (%.7f,%.7f)\n",points[p],points[partLength + p]);
            // }
            // for (int p=0; p<partLength; p++){
            //     printf("Distance from master before distribute: %.7f\n",distances[p]);
            // }
            MPI_Bcast(gtMedian,sub_size,MPI_INT,0,sub_comm);
            MPI_Bcast(ltMedian,sub_size,MPI_INT,0,sub_comm);
            for (int p=0; p<sub_size; p++){
                printf("Greater[%d]-%d from master before distribute: %d\n",p,k,gtMedian[p]);
                printf("Lower[%d]-%d from master before distribute: %d\n",p,k,ltMedian[p]);
            }
            printf("Distribution is starting\n");
            Distribute(points,ltMedian,gtMedian,sub_id,sub_size,partLength,sub_comm);       // Anakatanemoume ta points gia na plhrountai ta criteria mas
            printf("Distribution has been completed\n");
            calcDistance(points, partLength, vp, distances);                                // Ypologizoume ta nea distances
            checkSums(distances,points,partLength,median,gtMedian,eqMedian,ltMedian,sub_comm);
            // for (int p=0; p<partLength; p++){
            //     printf("Point from master after distribute: (%.7f,%.7f)\n",points[p],points[partLength + p]);
            // }
            // for (int p=0; p<partLength; p++){
            //     printf("Distance from master after distribute: %.7f\n",distances[p]);
            // }
            printf("\n\n");
            for (int p=0; p<sub_size; p++){
                printf("Greater[%d]-%d from master after distribute: %d\n",p,k,gtMedian[p]);
                printf("Lower[%d]-%d from master after distribute: %d\n",p,k,ltMedian[p]);
            }
            distrValidation(sub_id,sub_size,distances,median,partLength,sub_comm);      // Validation
            free(vp);
            free(distances);
            //free(distances_bk);
            MPI_Comm_free(&masters_comm);
        }
        else
        {
            float* vp = (float*)malloc(dimension*sizeof(float));
            MPI_Bcast(vp,dimension,MPI_FLOAT,0,sub_comm);
            float* distances = (float*)malloc(partLength*sizeof(float));
            //float* distances_bk = (float*)malloc(partLength*sizeof(float));
            calcDistance(points, partLength, vp, distances);
            //memcpy(distances_bk,distances,partLength*sizeof(float));
            // for (int p=0; p<partLength*dimension; p++){
            //     printf("Point from slave: %.7f\n",points[p]);
            // }
            // for (int p=0; p<partLength; p++){
            //     printf("Distance: %.7f\n",distances[p]);
            // }
            // printf("Vantage point from slave%d: %f + %f\n",k, vp[0], vp[1]);
            slavePart(sub_id,partLength,distances,size/(1<<k),points,sub_comm);
            MPI_Bcast(&median,1,MPI_FLOAT,0,sub_comm);
            checkSums(distances,points,partLength,median,gtMedian,eqMedian,ltMedian,sub_comm);
            //swapPoints(points,distances,distances_bk,partLength,dimension);
            MPI_Bcast(gtMedian,sub_size,MPI_INT,0,sub_comm);
            MPI_Bcast(ltMedian,sub_size,MPI_INT,0,sub_comm);
            // for (int p=0; p<partLength*dimension; p++){
            //     printf("Point from slave before distribute: %f\n",points[p]);
            // }
            // for (int p=0; p<partLength; p++){
            //     printf("Distance from slave before distribute: %f\n",distances[p]);
            // }
            printf("Distribution from slave is starting\n");
            Distribute(points,ltMedian,gtMedian,sub_id,sub_size,partLength,sub_comm);
            printf("Distribution from slave has been completed\n");
            calcDistance(points, partLength, vp, distances);
            checkSums(distances,points,partLength,median,gtMedian,eqMedian,ltMedian,sub_comm);
            // for (int p=0; p<partLength*dimension; p++){
            //     printf("Point from slave after distribute: %f\n",points[p]);
            // }
            // for (int p=0; p<partLength; p++){
            //     printf("Distance from slave after distribute: %f\n",distances[p]);
            // }
            distrValidation(sub_id,sub_size,distances,median,partLength,sub_comm);
            free(vp);
            free(distances);
            //free(distances_bk);
        }
        MPI_Comm_free(&sub_comm);
        
    }

    serialVP(points,localVPs,localMedians,partLength,partLength,1,processId);       // Ftiaxnoume ta topika trees


    MPI_Bcast(sharedVPs,noProcesses*dimension,MPI_FLOAT,0,MPI_COMM_WORLD);          // Bcast to koino dentro stous pantes
    MPI_Bcast(sharedMedians,noProcesses,MPI_FLOAT,0,MPI_COMM_WORLD);

    
    MPI_Barrier(MPI_COMM_WORLD);
    if (processId == 0){
        printf("Construction of tree has ended!\n");
        end_time = MPI_Wtime();
        double Worktime = end_time - start_time;
        printf("Worktime with MPI_Wtime: %f seconds\n",Worktime);
        for (int i=1; i<(noProcesses-1)*dimension; i++){
           // printf("Vantage Points: %f\n", sharedVPs[i]);
        }
        for (int i=1; i<noProcesses; i++){
            //printf("Medians: %f\n", sharedMedians[i]);
        }
    }
    else if (processId == 1){
        for (int i=1; i<noProcesses; i++){
            printf("Medians from pid1: %f\n", sharedMedians[i]);
        }
        for (int i=1; i<noProcesses; i++){
            printf("Vantage Points from pid1: (%f,%f)\n", sharedVPs[i],sharedVPs[noProcesses-1+i]);
        }
        // for (int i=1; i<partLength/2; i++){
        //     printf("Local Medians from pid1: %f\n", localMedians[i]);
        // }
    }
    
    free(sharedVPs);
    free(sharedMedians);
    free(points);
    MPI_Finalize();
    return 0;
}


/*========================FIND MEDIAN FUNCTIONS====================================
 * ================================================================================
 * ================================================================================
*/


/****Partitions the Array into larger and smaller than the pivot values****/
// Allagi gia float
void partition (float *array, float* points, int partLength, int elements, float pivot, float **arraysmall, float **arraybig,float** pointsBig,float** pointsSmall, int *endsmall, int *endbig)
{
    int right=elements-1;
    int left=0;
    int pos;
    if(elements==1)
    {
        if(pivot>array[0])
        {
            *endsmall=1;  //One value in the small part
            *endbig=0;   //Zero on the big one
            *arraysmall=array;   //There is no big array therefore NULL value
            *arraybig=NULL;
            *pointsSmall=points;
            *pointsBig=NULL;
                
        }
        else if(pivot<=array[0])
        {
            *endsmall=0;    //The exact opposite of the above actions.
            *endbig=1;
            *arraysmall=NULL;
            *arraybig=array;
            *pointsSmall=NULL;
            *pointsBig=points;
            
        }
    }
    else if(elements>1)
    {
        while(left<right)
        {
            while(array[left]<pivot)
            {
                left++;
                if(left>=elements)
                {
                    break;
                }
            }
            while(array[right]>=pivot)
            {
                right--;
                if(right<0)
                {
                    break;
                }
            }
            if(left<right)
            {
                swap_values(array,points,partLength,left,right);
            }
        }
        pos=right;
        if(pos<0)                   //Arrange the arrays so that they are split into two smaller ones.
        {                               //One containing the small ones. And one the big ones.
            *arraysmall=NULL;           //However these arrays are virtual meaning that we only save the pointer values of the beging and end of the "real" one.
            *pointsSmall = NULL;        
        }                               
        else
        {
            *arraysmall=array;
            *pointsSmall = points;
        }
        *endsmall=pos+1;
        *arraybig=&array[pos+1];
        *endbig=elements-pos-1;
        *pointsBig=&points[pos+1];
    }
}


/***==============================================***/
/***==============================================***/
/***=============SERIAL SELECTION==============***/
/***==============================================***/
/***==============================================***/

// Allagi gia float
float selection(float *array,int number,float *points,int partLength)
{
    float *arraybig;
    float *arraysmall;
    int endsmall=0;
    int endbig=0;
    float *arraytobeused;
    int i;
    int counter=0;
    int k;
    float* pointsBig;
    float* pointsSmall;
    float* pointstobeused;
    float pivot;
    float median;
    k=(int)number/2+1;
    arraytobeused=array;
    pointstobeused=points;
    for(;;)
    {
        pivot=arraytobeused[rand() % number];
        partition(arraytobeused,pointstobeused,partLength,number,pivot,&arraysmall,&arraybig,&pointsBig,&pointsSmall,&endsmall,&endbig);
        if(endbig>k)
        {
            number=endbig;
            arraytobeused=arraybig;
            pointstobeused=pointsBig;
            for(i=0;i<endbig;i++)
            {
                if(pivot==arraybig[i])
                    counter++;
                else
                    break;
            }
            if(counter==endbig)
            {
                median=arraybig[0];
                break;
            }
            else
                counter=0;
            //end of count equals
        }
        else if(endbig<k)
        {
            number=endsmall;
            arraytobeused=arraysmall;
            pointstobeused=pointsSmall;
            k=k-endbig;
        }
        else
        {
            median=pivot;
            break;
        }
    }
    return median;
}
