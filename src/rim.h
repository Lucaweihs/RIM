typedef struct List {
  struct ListNode* sentinal;
  struct ListNode* current;
  int length;
} List;

typedef struct ListNode {
  struct ListNode* prev;
  struct ListNode* next;
  void* value;
} ListNode;

typedef struct RimTree {
  struct RimNode* root;
  struct List* refList;
  int* refArray;
  int* discMat;
} RimTree;

typedef struct RimNode {
  struct RimNode* left;
  struct RimNode* right;
  struct RimNode* parent;
  double theta;
  ListNode* element;
} RimNode;

List* createEmptyList(void);
bool atEnd(List*);
bool next(List*);
void restart(List*);
int length(List*);
ListNode* firstNode(List*);
ListNode* lastNode(List*);
ListNode* currentNode(List*);
void appendNode(List*, ListNode*);
void* firstValue(List*);
void* lastValue(List*);
void* currentValue(List*);
void appendValue(List*, void*);

RimTree* createRimTree(RimNode*,List*);
//RimTree* createRimTreeFromReference(List*);
void canonicalPermutation(RimNode*);
//void logLikelihood(List*);


RimNode* createEmptyRimNode();
void attachLeft(RimNode*, RimNode*);
void attachRight(RimNode*, RimNode*);
void attachRankingToLeaves(RimNode*, List*);
int attachRankingToLeavesHelper(RimNode*, List*);

int countPartitions(int*, int, int, int, int, int, int);
int* discrepancyMat(int*, int, int);
int* rankingFromList(List*);