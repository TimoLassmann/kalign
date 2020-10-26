#include "tldevel.h"
#include "queue.h"

struct kchaos_tree{
        int** t;
        int* active;
        int chaos;
        int numseq;
};

/* The two key functions */
static int treeify_samples(struct kchaos_tree** tree,int numseq, int chaos);
static void traverse_kctree( struct kchaos_tree* tree, int index);

/* Memory allocation */
static int alloc_kchaos_tree(struct kchaos_tree** tree, int numseq, int chaos);
static void free_kchaos_tree(struct kchaos_tree* t);


struct node{
        struct node* l;
        struct node* r;
        int id;
};

struct topo{
        struct node** list;
        int index;
        int num_alloc;
};

static struct topo*  generate_tree_topologies(int n);
static void print_topo( struct node* n,int pos);

int main(int argc, char *argv[])
{
        struct kchaos_tree* tree = NULL;

        int numseq = 10;


        RUN(treeify_samples(&tree, numseq, 12));

        traverse_kctree(tree, numseq*2-2);

        free_kchaos_tree(tree);


        struct topo* t = NULL;
        t = generate_tree_topologies(7);
        int i;

        for(i = 0; i < t->index;i++){
                LOG_MSG("Tree %d:", i);
                print_topo(t->list[i], 0);
        }
        return OK;
ERROR:
        return FAIL;
}

/* Number of labelled tree topologies */


/* # A very simple representation for Nodes. Leaves are anything which is not a Node. */
/* class Node(object): */
/*   def __init__(self, left, right): */
/*     self.left = left */
/*     self.right = right */

/*   def __repr__(self): */
/*     return '(%s %s)' % (self.left, self.right) */

/* # Given a tree and a label, yields every possible augmentation of the tree by */
/* # adding a new node with the label as a child "above" some existing Node or Leaf. */
/* def add_leaf(tree, label): */
/*   yield Node(label, tree) */
/*   if isinstance(tree, Node): */
/*     for left in add_leaf(tree.left, label): */
/*       yield Node(left, tree.right) */
/*     for right in add_leaf(tree.right, label): */
/*       yield Node(tree.left, right) */

/* # Given a list of labels, yield each rooted, unordered full binary tree with */
/* # the specified labels. */
/* def enum_unordered(labels): */
/*   if len(labels) == 1: */
/*     yield labels[0] */
/*   else: */
/*     for tree in enum_unordered(labels[1:]): */
/*       for new_tree in add_leaf(tree, labels[0]): */
/*         yield new_tree */

struct topo* build_forest(struct topo* t, int* lab, int total_n,int index, int d)
{
        /* if(d == 0){ */
        /*         /\* create  *\/ */
        /* } */
        /* int i; */

        /* if(n == 1){ */
        /*         fprintf(stdout,"%d\n", lab[0]); */
        /* }else{ */
        /*         for(i = 0; i < n;i++){ */
        /*                 res = enum_unordered(lab+1,n-1); */
        /*         } */


        /* } */
}

void print_topo( struct node* n,int pos)
{
        if(n == NULL){
                for(int i=0;i<pos;i++){
                        fprintf(stdout," ");
                }
                fprintf(stdout,"*\n");

        }
        if(n->r){
        print_topo(n->r, pos+1);
        }


        for(int i=0;i<pos;i++){
                fprintf(stdout," ");
        }
        fprintf(stdout,"%d\n",n->id);
        if(n->l){
                print_topo(n->l, pos+1);
        }
}

struct topo*  generate_tree_topologies(int n)
{

        struct topo* res = NULL;

        struct topo* left = NULL;
        struct topo* right = NULL;

        struct node* node = NULL;

        int i,j,c;

        MMALLOC(res, sizeof(struct topo));
        res->list = NULL;
        res->num_alloc = 16;
        res->index = 0;
        MMALLOC(res->list, sizeof(struct node*) * res->num_alloc);
        res->list[0] = NULL;

        if( n % 2 == 0){
                LOG_MSG("Aborting: n = %d",n);
                return res;
        }
        if(n == 1){
                LOG_MSG("N == 1");
                MMALLOC(node, sizeof(struct node));
                node->l = NULL;
                node->r = NULL;
                node->id = 0;
                res->list[res->index] = node;

                res->index++;
                return res;
        }
        LOG_MSG("I have %d", n);
        for (i = 1; i < n; i += 2) {
                LOG_MSG("Split %d  %d",i , n-i-1);
                left  = generate_tree_topologies(i);
                right = generate_tree_topologies(n - i -1);

                for(j = 0; j < left->index;j++){
                        for(c = 0; c < right->index;c++){
                                node = NULL;
                                MMALLOC(node, sizeof(struct node));
                                node->l = left->list[j];
                                node->r = right->list[c];
                                node->id = 1;
                                res->list[res->index] = node;
                                res->index++;
                                if(res->index == res->num_alloc){
                                        res->num_alloc = res->num_alloc + res->num_alloc / 2;
                                        MREALLOC(res->list, sizeof(struct node*) * res->num_alloc);
                                }
                        }
                }
                MFREE(left->list);
                MFREE(left);
                MFREE(right->list);
                MFREE(right);

        }


        LOG_MSG("res %d", res->index);
        return res;
ERROR:
        LOG_MSG("Somethign ");
        exit(-1);

    /*     public List<TreeNode> allPossibleFBT(int N) { */
    /*     List<TreeNode> result = new ArrayList<>(); */
    /*     if (N % 2 == 0) { */
    /*         return result; */
    /*     } */
    /*     if (N == 1) { */
    /*         result.add(new TreeNode(0)); */
    /*         return result; */
    /*     } */
    /*     for (int i = 1; i < N; i += 2) { */
    /*         List<TreeNode> lefts = allPossibleFBT(i); */
    /*         List<TreeNode> rights = allPossibleFBT(N - i - 1); */
    /*         for (TreeNode l : lefts) { */
    /*             for (TreeNode r : rights) { */
    /*                 TreeNode root = new TreeNode(0); */
    /*                 root.left = l; */
    /*                 root.right = r; */
    /*                 result.add(root); */
    /*             } */
    /*         } */
    /*     } */
    /*     return result; */
    /* } */
}


/* Function to traverse the chaos-ary tree.
   1) check if all child nodes are "active", meaning they
   are either leaves or have been processed already.
   2) if a child in not active call traverse tree recursively
   3) continue into all nodes are visited
 */
void traverse_kctree( struct kchaos_tree* tree, int index)
{
        int** t = NULL;
        int* a = NULL;
        int i;
        int c;
        a = tree->active;
        t = tree->t;
        for(i = 0; i < tree->chaos;i++){
                c = t[index][i];
                if(!a[c] && c!= -1){
                        traverse_kctree(tree, c);
                }
        }
        LOG_MSG("Processing: ");
        for(i = 0; i < tree->chaos;i++){
                fprintf(stdout,"%d ", t[index][i]);
        }
        fprintf(stdout,"\n");

        /* I have processed this node  */
        tree->active[index] = 1;
}


int treeify_samples(struct kchaos_tree** tree,int numseq, int chaos)
{
        struct kchaos_tree* t = NULL;

        int* tmp = NULL;
        int i,j,c;
        int work;
        int cur_node;
        int bin_index;


        MMALLOC(tmp, sizeof(int)* chaos * 2);

        RUN(alloc_kchaos_tree(&t, numseq, chaos));

        for(j = 0; j <  chaos*2;j++){
                tmp[j] = -1;

        }

        queue q = q_new();
        for(i = 0; i < numseq;i++){
                enqueue(q, i+1);
        }

        /* print_queue(q); */

        work = 1;
        j= 0;
        cur_node = 0;
        bin_index = numseq;
        while(work != numseq*2 -2 ){
                /* print_queue(q); */
                if (!dequeue(q, &work)) {
                        break;
                }
                tmp[j] = work-1;
                /* t[cur_node][j] = work-1; */
                if(j){
                        tmp[j-1+chaos] = bin_index;
                        /* t[cur_node][j-1+chaos] = bin_index; */
                        bin_index++;
                }
                j++;
                if(j == chaos){
                        /* fprintf(stdout,"Finished: "); */
                        for(j = 0; j <  chaos*2;j++){
                                t->t[bin_index-1][j] = tmp[j];
                                /* fprintf(stdout,"%d ", tmp[j]); */
                                tmp[j] = -1;

                        }
                        /* fprintf(stdout,"  bin: %d \n", bin_index-1); */
                        enqueue(q, bin_index);
                        cur_node++;
                        j = 0;

                }
        }
        if(j  > 1){
                /* We have a half completed node.  */
                /* fprintf(stdout,"Finished: %d ",j); */
                for(j = 0; j <  chaos*2;j++){
                        t->t[bin_index-1][j] = tmp[j];
                        /* fprintf(stdout,"%d ", tmp[j]); */
                }
                /* fprintf(stdout,"  bin: %d \n", bin_index); */
        }
        /* for(i = 0; i < numseq*2;i++){ */
        /*         fprintf(stdout,"%d:  ", i); */
        /*         for(j = 0; j <  chaos*2;j++){ */
        /*                 fprintf(stdout,"%d ",t->t[i][j]); */
        /*         } */
        /*         fprintf(stdout,"\n"); */
        /* } */
        /* print_queue(q); */

        free_queue(q);

        MFREE(tmp);
        *tree = t;
        return OK;
ERROR:
        free_kchaos_tree(t);
        return FAIL;
}


int alloc_kchaos_tree(struct kchaos_tree** tree, int numseq, int chaos)
{
        struct kchaos_tree* t = NULL;
        int i;
        int j;

        MMALLOC(t, sizeof(struct kchaos_tree));

        t->numseq = numseq;
        t->chaos = chaos;
        t->active = NULL;
        t->t = NULL;

        RUN(galloc(&t->t,numseq*2,chaos*2));
        RUN(galloc(&t->active,numseq*2));

        for(i = 0; i < numseq;i++){
                t->active[i] = 1;
        }
        for(i = numseq; i < numseq*2;i++){
                t->active[i] = 0;
        }
        for(i = 0; i < numseq*2;i++){
                for(j = 0; j <  chaos*2;j++){
                        t->t[i][j] = 0;
                }
        }
        *tree = t;
        return OK;
ERROR:
        free_kchaos_tree(t);
        return FAIL;
}

void free_kchaos_tree(struct kchaos_tree* t)
{
        if(t){
                if(t->t){
                        gfree(t->t);
                }
                if(t->active){
                        gfree(t->active);
                }
                MFREE(t);
        }
}
