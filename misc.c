#define FLAGLENGTH 5
#define MAXMSGLENGTH 100
#define MAXFILENAMELENGTH 100

//================================================================
int showxerror(const char* buff) {
  //Prints our meesage provided as argument and stops the program
  fprintf(stderr,"%s\n",buff);
  exit(0);
}
//================================================================
//================================================================
int showerror(const char* buff) {
  //Prints our message provided as argument
  fprintf(stderr,"%s\n",buff);
  return 1;
}
//================================================================
//================================================================
int look_for_flag(int *ieof, char flag[FLAGLENGTH], FILE *unit)
{
  //This function looks for a flag which is a * followed by a string of length 
  //FLAGLENGTH. There maybe a number of blanks between the * and the string. 
  //They will be ignored. The returned flag is not string not including the *
  char last_char;
  int line_begins=1;
  int i;
  
  //As long as end of file has not been reached
  while(!*ieof){
    //we test for end of file
    *ieof=feof(unit);
    if(*ieof) return 0;  
    //reads one character
    fscanf(unit,"%c",&last_char);
    //skip any blank
    if (last_char==' ') continue;
    //At this point we found a character that is not a '' 
    if (line_begins){
      line_begins=0;
      if (last_char=='*') {
	last_char=' ';
	//skips all space following the '*' 
	while(last_char==' ') {
	  *ieof=feof(unit);
	  if(*ieof) return 0;  
	  fscanf(unit,"%c",&last_char);
	}
	//We found the first character of a flag
	flag[0]=last_char;
	//and we read the others
	for (i=1;i<FLAGLENGTH;i++) {
	  *ieof=feof(unit);
	  if(*ieof) return 0;  
	  fscanf(unit,"%c",&flag[i]);
	}
	/*we found a flag and now return*/
	return 0;
      }
    } 
    /*if we found a carriage return we are going to start a new line*/
    if (last_char=='\n') line_begins=1;
  }
  return 0;
}
//================================================================
//================================================================
int get_word(FILE *funit, char name[MAXFILENAMELENGTH], int *ieof)
//Read a word defined a a string with no blanks. The first encountered 
//blank will indicate the end of the string. 
{
  char last_char;          /*Last character read on the line*/ 
  int fnamlen=0;           /*File name length*/
  
  while(!*ieof){
    /*we test for end of file*/
    *ieof=feof(funit);
    if(*ieof) return 0;  
    /*reads one character*/
    fscanf(funit,"%c",&last_char);
    
    /*skip any blank before the name*/
    if (last_char==' ' && fnamlen==0) continue;
    /*any blank encountered when reading the name indicates 
      the end of the name*/
    if (last_char==' ' || last_char=='\n') break;
    
    name[fnamlen]=last_char;
    fnamlen=fnamlen+1;
    if (fnamlen>MAXFILENAMELENGTH) 
      return showerror("Word too long");
  }
  name[fnamlen]='\0';
  return 0;   
}
//================================================================
