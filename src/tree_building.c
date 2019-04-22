#include "tree_building.h"

#define NODESIZE 16

/* small hash implementation */
struct bignode{
        struct bignode *next;
        unsigned int pos[NODESIZE];
        unsigned int num;
};

struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos);
void big_remove_nodes(struct bignode *n);
void big_print_nodes(struct bignode *n);

/* distance calculation */
float** dna_distance(struct alignment* aln,struct parameters* param, int nj);
float dna_distance_calculation(struct bignode* hash[],int* p,int seqlen,int diagonals,float mode);


float** protein_wu_distance(struct alignment* aln,struct parameters* param, int nj);
float protein_wu_distance_calculation(struct bignode* hash[],const int* seq,const int seqlen,const int diagonals,const float mode);

/* tree building stuff */

struct aln_tree_node{
        struct aln_tree_node** links;
        int* internal_lables;
        int* path;
        int* profile;
        int* seq;
        int len;
        int done;
        int num;
};

struct aln_tree_node* real_upgma(float **dm,int ntree,int numseq);
struct aln_tree_node* real_nj(float** dm,int ntree, int numseq);

int* readtree(struct aln_tree_node* p,int* tree);
void ntreeify(struct aln_tree_node* p,int ntree);


int build_tree(struct alignment* aln,struct parameters* param, struct aln_param* ap)
{
        struct aln_tree_node* tree2 = NULL;
        float** dm = NULL;
        int* tree = NULL;
        int numseq;
        int a,c;
        int i;

        ASSERT(aln != NULL, "No alignment.");
        ASSERT(param != NULL, "No input parameters.");
        ASSERT(ap != NULL, "No alignment parameters.");

        tree = ap->tree;
        numseq = aln->numseq;

        /* if we don't build a tree fill and return  */
        if(param->ntree <= 1){
                tree[0] = 0;
                tree[1] = 1;

                c = numseq;
                tree[2] = c;
                a = 2;
                for ( i = 3; i < (numseq-1)*3;i+=3){
                        tree[i] = c;
                        tree[i+1] = a;
                        c++;
                        tree[i+2] = c;
                        a++;
                }
                return OK;
        }

        /* calculate distances  */



        if(param->dna == 1){
                if(byg_start(param->tree,"njNJ") != -1){
                        RUNP(dm =  dna_distance(aln,param,1));
                }else{
                        RUNP(dm =  dna_distance(aln,param,0));
                }
        }else{
                if(byg_start(param->tree,"njNJ") != -1){
                        RUNP(dm =  protein_wu_distance(aln,param,1));
                }else{
                        RUNP(dm =  protein_wu_distance(aln,param,0));
                }
        }

        if(byg_start(param->tree,"njNJ") != -1){
                tree2 = real_nj(dm,param->ntree,aln->numseq);
        }else{
                tree2 = real_upgma(dm,param->ntree,aln->numseq);
        }

        gfree(dm);

        //if(param->print_tree){
        //print_tree(tree2,aln,param->print_tree);
        //}


        tree[0] = 1;

        if(param->ntree > 2){
                ntreeify(tree2,param->ntree);
        }else{
                tree = readtree(tree2,tree);
                for (i = 0; i < (numseq*3);i++){
                        tree[i] = tree[i+1];
                }
                free(tree2->links);
                free(tree2->internal_lables);
                free(tree2);
        }



        return OK;
ERROR:
        return FAIL;
}




float** protein_wu_distance(struct alignment* aln,struct parameters* param, int nj)
{
        struct bignode* hash[1024];
        float** dm = NULL;
        int*p =0;
        int i,j,a,b;
        unsigned int hv;
        int numseq;
        int numprofiles;
        float min;
        float cutoff;

        ASSERT(aln != NULL,"No alignment");

        numseq = aln->numseq;
        numprofiles = aln->num_profiles;

        for (i = 0;i < 1024;i++){
                hash[i] = 0;
        }

        i = numseq;
        if (nj){
                i = numprofiles;
        }

        RUNP(dm = galloc(dm,i,i,0.0f));
        fprintf(stderr,"Distance Calculation:\n");
        b = (numseq*(numseq-1))/2;
        a = 1;

        for (i = 0; i < numseq-1;i++){
                p = aln->s[i];

                for (j = aln->sl[i]-2;j--;){
                        //for(j = 0; j < si->sl[i]-2;j++){
                        //hv = (p[j+1] << 5) + p[j+2];
                        //hash[hv] = big_insert_hash(hash[hv],j);
                        hv = (p[j] << 5) + p[j+1];
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = (p[j] << 5) + p[j+2];
                        hash[hv] = big_insert_hash(hash[hv],j);
                }
                for (j = i+1; j < numseq;j++){
                        min =  (aln->sl[i] > aln->sl[j]) ? aln->sl[j] :aln->sl[i];
                        cutoff = param->internal_gap_weight *min + param->zlevel;
                        //cutoff = param->zlevel;
                        p = aln->s[j];
                        dm[i][j] = protein_wu_distance_calculation(hash,p,aln->sl[j],aln->sl[j]+aln->sl[i],cutoff);
                        //fprintf(stderr,"%d-%d:%f\n",i,j,dm[i][j]);
                        //exit(0);
                        //dm[i][j] /= min;
                        //dm[i][j] /= (si->sl[i] > si->sl[j]) ? si->sl[j] :si->sl[i];
                        dm[j][i] = dm[i][j];
                        fprintf(stderr,"\r%8.0f percent done",(float)a /(float)b * 100);
                        a++;
                }


                for (j = 1024;j--;){
                        if (hash[j]){
                                big_remove_nodes(hash[j]);
                                hash[j] = 0;
                        }
                }
        }
        return dm;
ERROR:
        return NULL;
}


float protein_wu_distance_calculation(struct bignode* hash[],const int* seq,const int seqlen,const int diagonals,const float mode)
{

        struct bignode* node_p;
        unsigned int* d = NULL;
        unsigned int* tmp = NULL;
        float out = 0.0;
        register int i,j;
        register int c;
        register int num;
        register unsigned int hv;

        d = malloc(sizeof(unsigned int)*diagonals);
        //for (i = diagonals;i--;){
        for (i = 0;i < diagonals;i++){
                d[i] = 0;
        }
        for (i = seqlen-2;i--;){
                //for(i = 0; i < seqlen-2;i++){
                /*hv = (seq[i+1] << 5) + seq[i+2];

                  node_p = hash[hv];
                  while(node_p){
                  tmp = node_p->pos;
                  for(j = 0;j < node_p->num;j++){
                  d[tmp[j]]++;
                  }
                  node_p = node_p->next;
                  }*/
                hv = (seq[i] << 5) + seq[i+1];
                //printf("3:%d\n",hv);
                node_p = hash[hv];
                while(node_p){
                        tmp = node_p->pos;
                        num = node_p->num;
                        for(j = 0;j < num;j++){
                                c = tmp[j];
                                d[c]++;
                                c++;
                                d[c]++;
                        }
                        node_p = node_p->next;
                }
                hv = (seq[i] << 5) + seq[i+2];

                node_p = hash[hv];

                while(node_p){
                        tmp = node_p->pos;
                        num = node_p->num;
                        for(j = 0;j < num;j++){
                                c = tmp[j];
                                d[c]++;
                        }
                        node_p = node_p->next;
                }
                d++;
        }
        //exit(0);
        d -= (seqlen-2);

        for (i = diagonals;i--;){
                //d[i] /= minlen;

                //fprintf(stderr,"%d ",d[i]);
                if(d[i] > mode){
                        out += d[i];
                        //	printf("%f	%d\n",d[i]/ minlen,d[i]);
                }
        }
        free(d);
        return out;
}



float** dna_distance(struct alignment* aln,struct parameters* param, int nj)
{
        struct bignode* hash[1024];
        float** dm = NULL;
        int *p = 0;
        int i,j,a,b;
        unsigned int hv;
        int numseq;
        int numprofiles;

        ASSERT(aln != NULL,"No alignment");

        numseq = aln->numseq;
        numprofiles = aln->num_profiles;

        fprintf(stderr,"Distance Calculation:\n");


        for (i = 0;i < 1024;i++){
                hash[i] = 0;
        }

        i = numseq;
        if (nj){
                i = numprofiles;
        }

        RUNP(dm = galloc(dm,i,i,0.0f));


        b = (numseq*(numseq-1))/2;
        a = 1;

        for (i = 0; i < numseq-1;i++){
                p = aln->s[i];
                for (j = aln->sl[i]-5;j--;){
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+3]&3)<<2) + (p[j+4]&3);//ABCDE
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+3]&3)<<2) + (p[j+5]&3);//ABCDF
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+4]&3)<<2) + (p[j+5]&3);//ABCEF
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+3]&3)<<4)  + ((p[j+4]&3)<<2) + (p[j+5]&3);//ABDEF
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+2]&3)<<6) + ((p[j+3]&3)<<4) + ((p[j+4]&3)<<2) + (p[j+5]&3);//ACDEF
                        hash[hv] = big_insert_hash(hash[hv],j);
                }
                for (j = i+1; j < numseq;j++){


                        //min =  (aln->sl[i] > aln->sl[j]) ?aln->sl[j] :aln->sl[i];
                        dm[i][j] = dna_distance_calculation(hash,aln->s[j],aln->sl[j],aln->sl[j]+aln->sl[i],param->zlevel);
                        dm[i][j] /= (aln->sl[i] > aln->sl[j]) ?aln->sl[j] :aln->sl[i];
                        dm[j][i] = dm[i][j];
                        fprintf(stderr,"\r%8.0f percent done",(float)a /(float)b * 100);
                        a++;
                }

                for (j = 1024;j--;){
                        if (hash[j]){
                                big_remove_nodes(hash[j]);
                                hash[j] = 0;
                        }
                }
        }
        return dm;
ERROR:
        return NULL;
}

float dna_distance_calculation(struct bignode* hash[],int* p,int seqlen,int diagonals,float mode)
{

        struct bignode* node_p;
        float out = 0.0;
        unsigned int* tmp = NULL;
        unsigned int* d = NULL;
        int i,j;
        unsigned int hv;

        d = malloc(sizeof(int)*diagonals);
        for (i = 0;i < diagonals;i++){
                d[i] = 0;
        }
        for (i = seqlen-5;i--;){

                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+3]&3)<<2) + (p[i+4]&3);//ABCDE
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }


                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+3]&3)<<2) + (p[i+5]&3);//ABCDF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }
                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+4]&3)<<2) + (p[i+5]&3);//ABCEF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }
                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+3]&3)<<4)  + ((p[i+4]&3)<<2) + (p[i+5]&3);//ABDEF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }
                hv = ((p[i]&3)<<8) + ((p[i+2]&3)<<6) + ((p[i+3]&3)<<4) + ((p[i+4]&3)<<2) + (p[i+5]&3);//ACDEF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }

                d++;
        }
        //exit(0);
        d -= (seqlen-5);

        for (i = diagonals;i--;){
                //d[i] /= minlen;

                //printf("%d ",d[i]);

                if(d[i] > mode){
                        //fprintf(stderr,"%f	%d\n",d[i]/ minlen,d[i]);
                        out += d[i];
                }
        }
        free(d);
        return out;
}



struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos)
{
        struct bignode* p = NULL;
        if(n){
                if(n->num < NODESIZE){
                        n->pos[n->num] = pos;
                        n->num++;
                        return n;
                }else{
                        MMALLOC(p, sizeof(struct bignode));
                        p->pos[0] = pos;
                        p->num = 1;
                        p->next = n;
                }
        }else{
                MMALLOC(p, sizeof(struct bignode));
                p->pos[0] = pos;
                p->num = 1;
                p->next = n;
        }
        return p;
ERROR:
        return NULL;
}

void big_remove_nodes(struct bignode *n)
{
        struct bignode* p = NULL;
        while(n){
                p = n;
                n = n->next;
                MFREE(p);
        }
}

void big_print_nodes(struct bignode *n)
{
        int i;
        while(n){
                for (i = 0; i < n->num;i++){
                        fprintf(stderr,"%d ",n->pos[i]);
                }
                n = n->next;
        }
}



struct aln_tree_node* real_nj(float** dm,int ntree, int numseq)
{
        int i,j;
        //float **dm = 0;
        float *r = NULL;
        float *r_div = NULL;
        int *active = NULL;
        int node = 0;
        float min = 0;
        int join_a = 0;
        int join_b = 0;
        int leaves = 0;

        struct aln_tree_node** tree = NULL;
        struct aln_tree_node* tmp = NULL;

        leaves = numseq;


        MMALLOC(r,(numseq*2-1) *sizeof(float));
        MMALLOC(r_div,(numseq*2-1) *sizeof(float));
        MMALLOC(active,(numseq*2-1)*sizeof(int));

        for ( i = 0;i < numseq*2-1;i++){
                active[i] = 0;
        }
        for ( i = 0;i < numseq;i++){
                active[i] = 1;
        }


        MMALLOC(tree,sizeof(struct aln_tree_node*)*(numseq*2-1));
        for (i=0;i < numseq*2-1;i++){
                tree[i] = NULL;
                MMALLOC(tree[i],sizeof(struct aln_tree_node));
                tree[i]->done = 1;
                tree[i]->num = i;
                tree[i]->path = 0;
                tree[i]->profile = 0;
                tree[i]->seq = 0;//seq[i];
                tree[i]->len = 0;//len[i];
                tree[i]->internal_lables = NULL;

                MMALLOC(tree[i]->internal_lables,sizeof(int)*(ntree+(ntree-1)));
                MMALLOC(tree[i]->links,sizeof(struct aln_tree_node*)*(ntree+(ntree-1)));

                for ( j =0;j < (ntree+(ntree-1));j++){
                        tree[i]->links[j] = 0;
                        tree[i]->internal_lables[j] = 0;
                }
        }

        node = numseq;
        while (node != numseq*2 -1){
                for (i = 0;i<numseq*2-1;i++){
                        if (active[i]){
                                r[i] = 0;
                                for (j = 0;j < numseq*2-1;j++){
                                        if (active[j]){
                                                r[i] += (i<j) ?dm[i][j]:dm[j][i];
                                        }
                                }
                                r_div[i] = r[i] / (leaves-2);
                        }
                }
                for ( j = 0;j < numseq*2-1;j++){
                        if (active[j]){
                                for ( i = j+1;i < numseq*2-1;i++){
                                        if (active[i]){
                                                dm[i][j] = dm[j][i] - (r[i] + r[j])/2;
                                        }
                                }
                        }
                }
                min = -FLT_MAX;
                for ( j = 0;j < numseq*2-1;j++){
                        if (active[j]){
                                for ( i = j+1;i < numseq*2-1;i++){
                                        if (active[i]){
                                                if (dm[i][j] > min){
                                                        min = dm[i][j];
                                                        join_a = j;
                                                        join_b = i;
                                                }
                                        }
                                }
                        }
                }
                //join_a always smaller than join_b && both smaller than node
                dm[join_a][node] =  dm[join_a][join_b]/2 + (r_div[join_a] - r_div[join_b])/2;
                dm[join_b][node] =  dm[join_a][join_b] - dm[join_a][node];

                tree[node]->num = node;
                tree[node]->links[0] = tree[join_a];
                tree[node]->links[1] = tree[join_b];
                tree[node]->internal_lables[0] = node;
                tree[node]->internal_lables[1] = 0;


                active[join_a] = 0;
                active[join_b] = 0;

                for (i = 0;i<numseq*2-1;i++){
                        if (active[i]){
                                dm[i][node] = (i>join_a) ? dm[join_a][i]: dm[i][join_a];
                                dm[i][node] -= dm[join_a][node];
                                dm[i][node] += (i > join_b) ? dm[join_b][i] : dm[i][join_b] ;
                                dm[i][node] -= dm[join_b][node];
                                dm[i][node] /= 2;
                        }
                }
                active[node] = 1;
                node++;
        }


        MFREE(r);
        MFREE(r_div);
        MFREE(active);
        tmp = tree[node-1];
        MFREE(tree);
        return tmp;
ERROR:
        return NULL;
}

struct aln_tree_node* real_upgma(float **dm,int ntree,int numseq)
{
        int i,j;
        int *as = NULL;
        float max;
        int node_a = 0;
        int node_b = 0;
        int cnode = numseq;
        int numprofiles;

        struct aln_tree_node** tree = NULL;
        struct aln_tree_node* tmp =NULL;

        numprofiles = (numseq << 1) - 1;

        MMALLOC(as,sizeof(int)*numseq);
        for (i = numseq; i--;){
                as[i] = i+1;
        }

        MMALLOC(tree,sizeof(struct aln_tree_node*)*numseq);
        for (i=0;i < numseq;i++){
                tree[i] = NULL;
                MMALLOC(tree[i],sizeof(struct aln_tree_node));
                tree[i]->done = 1;
                tree[i]->num = i;
                tree[i]->path = 0;
                tree[i]->profile = 0;
                tree[i]->seq = 0;//seq[i];
                tree[i]->len = 0;//len[i];
                /*
                  Needs to be +2 because:
                  at n = 3 is is possible to get a perfectly balanced binary tree with 4 sequences at intermediate nodes
                */
                /*tree[i]->links = malloc(sizeof(struct aln_tree_node*)*2);

                  for ( j =0;j < 2;j++){
                  tree[i]->links[j] = 0;
                  }*/
                tree[i]->internal_lables = NULL;
                tree[i]->links = NULL;
                MMALLOC(tree[i]->internal_lables,sizeof(int)*(ntree+(ntree-1)));
                MMALLOC(tree[i]->links,sizeof(struct aln_tree_node*)*(ntree+(ntree-1)));

                for ( j =0;j < (ntree+(ntree-1));j++){
                        tree[i]->links[j] = 0;
                        tree[i]->internal_lables[j] = 0;
                }
        }

        while (cnode != numprofiles){
                max = -FLT_MAX;
                for (i = 0;i < numseq-1; i++){
                        if (as[i]){
                                for ( j = i + 1;j < numseq;j++){
                                        if (as[j]){
                                                if (dm[i][j] > max){
                                                        max = dm[i][j];
                                                        node_a = i;
                                                        node_b = j;
                                                }
                                        }
                                }
                        }
                }
                tmp = NULL;
                MMALLOC(tmp,sizeof(struct aln_tree_node));
                tmp->done = 0;
                tmp->path = 0;
                tmp->profile = 0;
                tmp->num = cnode;
                tmp->seq = 0;
                tmp->len = 0;
                tmp->links =NULL;
                tmp->internal_lables = NULL;
                MMALLOC(tmp->links,sizeof(struct aln_tree_node*)*(ntree+(ntree-1)));
                MMALLOC(tmp->internal_lables,sizeof(int)*(ntree+(ntree-1)));


                tmp->links[0] = tree[node_a];
                tmp->links[1] = tree[node_b];
                tmp->internal_lables[0] = cnode;
                tmp->internal_lables[1] = 0;

                for ( i =2;i < (ntree+(ntree-1));i++){
                        tmp->links[i] = 0;
                        tmp->internal_lables[i] = 0;

                }


                tree[node_a] = tmp;
                tree[node_b] = 0;

                /*deactivate  sequences to be joined*/
                as[node_a] = cnode+1;
                as[node_b] = 0;
                cnode++;

                /*calculate new distances*/
                for (j = numseq;j--;){
                        if (j != node_b){
                                dm[node_a][j] = (dm[node_a][j] + dm[node_b][j])*0.5;
                        }
                }
                dm[node_a][node_a] = 0.0f;
                for (j = numseq;j--;){
                        dm[j][node_a] = dm[node_a][j];
                        dm[j][node_b] = 0.0f;
                        dm[node_b][j] = 0.0f;
                }
        }
        tmp = tree[node_a];



        MFREE(tree);
        MFREE(as);
        return tmp;
ERROR:
        return NULL;
}

int* readtree(struct aln_tree_node* p,int* tree)
{
        if(p->links[0]){
                tree = readtree(p->links[0],tree);
        }
        if(p->links[1]){
                tree = readtree(p->links[1],tree);
        }

        if(p->links[0]){
                if(p->links[1]){
                        tree[tree[0]] = p->links[0]->num;
                        tree[tree[0]+1] = p->links[1]->num;
                        tree[tree[0]+2] = p->num;
                        tree[0] +=3;
                        MFREE(p->links[0]->internal_lables);
                        MFREE(p->links[0]->links);
                        MFREE(p->links[0]);
                        MFREE(p->links[1]->internal_lables);
                        MFREE(p->links[1]->links);
                        MFREE(p->links[1]);
                }
        }
        return tree;
}


void ntreeify(struct aln_tree_node* p,int ntree)
{
        int i = 0;
        int c = 0;
        struct aln_tree_node* tmp1 = 0;
        struct aln_tree_node* tmp2 = 0;
        if (p->links[0]){
                ntreeify(p->links[0],ntree);
        }
        if (p->links[1]){
                ntreeify(p->links[1],ntree);
        }

        if (!p->done){
                tmp1 = p->links[0];
                tmp2 = p->links[1];

                p->done = tmp1->done + tmp2->done;
                i = 0;
                c = 0;
                if(tmp1->done != 1){

                        while(tmp1->internal_lables[i]){
                                p->internal_lables[c] = tmp1->internal_lables[i];
                                i++;
                                c++;
                        }
                        if(tmp2->done != 1){
                                i = 0;
                                while(tmp2->internal_lables[i]){
                                        p->internal_lables[c] = tmp2->internal_lables[i];
                                        c++;
                                        i++;
                                }
                        }
                }else if(tmp2->done != 1){
                        i = 0;
                        while(tmp2->internal_lables[i]){
                                p->internal_lables[c] = tmp2->internal_lables[i];
                                c++;
                                i++;
                        }
                }
                p->internal_lables[c] = p->num;

                //fprintf(stderr,"%d:%d	%d:%d		%d\n",tmp1->num,tmp1->internal_lables[0],tmp2->num,tmp2->internal_lables[0],p->num);
                /*for (i = 0; i< c;i++){
                  fprintf(stderr,"il:%d ",p->internal_lables[i]);
                  }
                  fprintf(stderr,"\n");*/


                if (tmp1->done > 1){
                        for ( i = 0;i < tmp1->done;i++){
                                p->links[i] = tmp1->links[i];
                                tmp1->links[i] = 0;
                        }
                }

                if (tmp2->done > 1){
                        for ( i = 0; i < tmp2->done;i++){
                                p->links[tmp1->done+i] = tmp2->links[i];
                                tmp2->links[i] = 0;
                        }
                        MFREE(tmp2->internal_lables);
                        MFREE(tmp2->links);
                        MFREE(tmp2);
                }else{
                        p->links[tmp1->done] = tmp2;
                }
                //	fprintf(stderr,"p->num:%d\n",p->num);
                p->links[p->done] = 0;

                if (tmp1->done > 1){
                        MFREE(tmp1->internal_lables);
                        MFREE(tmp1->links);
                        MFREE(tmp1);
                }

                if (p->done >= ntree){
                        p->done = 1;
                        /*i = 0;
                          while(p->internal_lables[i]){
                          i++;
                          }
                          p->internal_lables[i] = p->num;*/
                }
        }
}
