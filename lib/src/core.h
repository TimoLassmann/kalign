#ifndef CORE_H
#define CORE_H

#ifdef CORE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


/* idea:
   use a shared struct holding dm & msa & empty root

   Process:
   put root on task list

   if (node is ready to align - do so )

   }else{


   while # tries < x
   take task, run split , update L / R arrays if better than previous split

   then add left and right node

   contine.
   }

   if(nseq at node == 1){
   send possible align signal up;
   if(align counter ==2 ){
   add parent node to queue

   }

   }
   if (number of left == 1 && number of right ==1){
   send possible align to

   }



*/


#undef CORE_IMPORT
#undef EXTERN


#endif
