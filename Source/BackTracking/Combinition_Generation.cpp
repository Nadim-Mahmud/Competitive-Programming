/**
*   Name : All combination generation
*   This code will print all nCr
*   While ara is elements of n
*/

int ara[20],sln[20],n,r;

void combinition(int pos,int last){
    int i;
    if(pos>r){
        for(i=1;i<=r;i++){
            if(i==r) printf("%d\n",sln[i]);
            else printf("%d ",sln[i]);
        }
        return;
    }

    for(i=last+1;i<=n;i++){
        sln[pos] = ara[i];
        combinition(pos+1,i);
    }
}

