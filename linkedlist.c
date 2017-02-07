

struct lnklstdata {
  float time;  
  float ampl;
};

struct node {
  struct lnklstdata data;
  struct node *prev;
  struct node *next;
};

struct list {
  struct node *head;
  struct node *tail;
  struct node *current;
  int nbnodes;
};
int initialize_list(struct list *l);
//=============================================================================
int initialize_list(struct list *l){
  l->nbnodes=0;
  l->head=NULL;
  l->tail=NULL;
  l->current=NULL;
  return 0;
}
//=============================================================================
//=============================================================================
int makenode(struct list *l){
  struct node *newnode;

  //Allocate memory for the new node
  if ((newnode  = (struct node *)malloc(sizeof(struct node)))==NULL) {
    fprintf(stderr,"\n\n Memory unavailable for more nodes\n");
    fprintf(stderr,"Exiting from Code :");
    exit(0);
  }
  //Are we creating the first node or are we adding a node? 
  if(l->current==NULL) {
    l->head=newnode;
    l->tail=newnode;
    l->current=newnode;
    newnode->prev=NULL;
    newnode->next=NULL;
  } else {
    newnode->prev=l->tail;
    newnode->next=NULL;    
    newnode->prev->next=newnode;
    l->tail=newnode;
    l->current=newnode;
  }
  l->nbnodes++;
  return 0; 
}
//=============================================================================
//=============================================================================
int reset_current(struct list *l){
  l->current = l-> head;
  return 0;
}
//=============================================================================
//=============================================================================
int next_one(struct list *l){
  if(l->current==NULL) return -1;
  if(l->current->next ==NULL) return -1;
  else l->current = l->current->next;
  return 0;
}
//=============================================================================
//=============================================================================
int delete_current(struct list *l){

  struct node *tmp;
  tmp=l->current;

  //Make sure there is a node to delete
  if (tmp == NULL) return -1;

  //If current node is the head and the tail
  if(tmp->prev == NULL && tmp->next == NULL) {
    l->head =NULL;
    l->tail=NULL;
    l->current=NULL;
    free(tmp);
    l->nbnodes--; 
    return 0;
  }

  //If current node is the list head
  if(tmp->prev == NULL) {
    l->head = tmp->next;
    l->current=tmp->next;
    l->head->prev=NULL;
    free(tmp);
    l->nbnodes--; 
    return 0;
  }

  //If current node is the list tail
  if(tmp->next == NULL) {
    l->tail = tmp->prev;
    l->current=tmp->prev;
    l->tail->next=NULL;
    free(tmp);
    l->nbnodes--;
    return 0;
  }

  //If we get here, the curent node has one before and one after
  tmp->next->prev=tmp->prev;
  tmp->prev->next=tmp->next;
  l->current=tmp->next;
  free(tmp);
  l->nbnodes--;
  return 0;
}
//=============================================================================
//=============================================================================
int clearlist(struct list *l){
  if(l->nbnodes==0) return 0;
  reset_current(l);
  do{
    delete_current(l);
  }  while(l->nbnodes>0);
  
  return 0;
}
//=============================================================================
//=============================================================================
