#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <queue>
#include <map>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include "RTree.h"

/// dUnit KNN查询半径扩张速度
#define dUnit 0.0001
using namespace std;
double b_lo, b_la;
int qid;


/// 定义了一个用于插入RTREE的节点，为矩形；如果要求插入点的话，可以把矩形的长宽设为0
template<class ElemType>
struct Rect
{
  Rect()  {}

  Rect(ElemType a_minX, ElemType a_minY, ElemType a_maxX, ElemType a_maxY)
  {
    min[0] = a_minX;
    min[1] = a_minY;

    max[0] = a_maxX;
    max[1] = a_maxY;
  }


  ElemType min[2];
  ElemType max[2];
};

///这是RTREE中的Node节点，为RTREE的DATATYPE，每个RTREE的节点均拥有一个Node作为指示标记
/// lo表示坐标的经度 la表示坐标的纬度
/// t_lo、t_la 表示在进行KNN时，该点被哪个query点包括了 t_lo为queryNode的经度，t_la为queryNode的纬度
/// query_index表示queryNode在queryPoints中的下标
/// route_index表示该点属于的路径的集合，如果route_index.empty(),那么该点不属于任何路径
struct Node
{
  Node() {}

  Node(double a_lo, double a_la){
    lo = a_lo;
    la = a_la;
  }
  vector<int> route_index;
  int query_index = 0;
  double lo = 0;
  double la = 0;
  double t_lo = 0;
  double t_la = 0;
};

vector<Node> tmpTargetNode;

/// 读取文件中的经纬度
/// 注意到：经纬度每相差一度，那么在实际地理位置上大约相差111km
/// 对于长距离的路径匹配（城市与城市之间），IKNN可能并不是一个高效的匹配算法
/// 所以，在本次实验中，为了控制地理位置的差距
/// 选择经度整数部分为-118，纬度整数部分为33的点
vector<Node> LoadData(char* filename, int& num){
    fstream f;
    double point_lo, point_la;
    char tmpchar;
    vector<Node> v;
    f.open(filename, ios::in);
    if(f){
        cout << "Loading......" << endl;
    }
    while(!f.eof()){
        f >> point_lo;
        f >> tmpchar;
        f >> point_la;
        if(floor(point_lo) == 33 && floor(point_la) == -119){
          Node tmp(point_lo, point_la);
          v.push_back(tmp);
          num++;
        }
    }
    return v;
}

/// 读取进来的点构造RTREE， 每个RTREE的DataType为Node类型，方便IKNN算法的实现
void BuildRtree(vector<Node> points, RTree<Node, double, 2, double>& myTree){
  for (int i = 0; i < points.size(); i++){
    Rect<double> *tmp = new Rect<double>(points[i].lo, points[i].la, points[i].lo, points[i].la);
    myTree.Insert(tmp->min, tmp->max, points[i]);
  }
}

/// 为RTREE 命中查询的回掉函数
/// 每当对一个queryPoint进行KNN时，如果命中，那么就将命中的点的DataType（Node）加入到vector中储存起来
/// 并修改该点的t_la，t_lo，query_index值
bool mySearchCallBack(Node node){
  node.t_la = b_la;
  node.t_lo = b_lo;
  node.query_index = qid;
  tmpTargetNode.push_back(node);
  return true;
}

/// 计算距离，为e的负欧氏距离
double calculateDist(Node a){
  return exp(-sqrt(pow(a.la - a.t_la, 2) + pow(a.lo - a.t_lo, 2)));
}

/// 用于sort的比价函数，选择一个点而不选择另外一个点
/// 是因为被选中的点距离queryPoint的距离更近
/// 如果距离更近，那么按照我们的算法计算出来的距离相似度就越大
bool myComp(Node a, Node b){
  return calculateDist(a) > calculateDist(b);
}

// 降序
bool myfunction(pair<int, double> m1, pair<int, double> m2){
  return m1.second > m2.second;
}

double calDistWithLoLa(double q_lo, double q_la, Node a){
  return exp(-sqrt(pow(a.lo - q_lo, 2) + pow(a.la - q_la, 2)));
}

double calculateSim(vector<vector<double>> queryPoints, vector<Node> path, map<int, double> simPath){
  double result = 0;
  for (int i = 0; i < queryPoints.size(); i++){ //对于每个queryPoint
    if(simPath.find(i) == simPath.end()){ //如果该点没有在路径遍历范围中
      double tmp = calDistWithLoLa(queryPoints[i][0], queryPoints[i][1], path[0]);
      for (int j = 1; j < path.size(); j++){
        if(calDistWithLoLa(queryPoints[i][0], queryPoints[i][1], path[j]) > tmp)
          tmp = calDistWithLoLa(queryPoints[i][0], queryPoints[i][1], path[j]);
      }
      result += tmp;
    }
    else{
      result += simPath[i];
    }
  }
  return result;
}

/// IKNN中的KNN算法
/// 由于一共需要命中k个点
/// 一个简单的算法就是 依次扩大搜索的矩阵长宽，直到命中k个点
/// 如果由于一次长宽的扩大导致命中了多于k个的点，那么取最近的k个点输出,即相似距离最大的点输出
vector<Node> qKNN(int queryIndex, int k, RTree<Node, double, 2, double> myTree, double q_lo, double q_la){
  int totalNum = 0;
  int num = 0;
  double iterTime = 1;
  b_lo = q_lo;
  b_la = q_la;
  qid = queryIndex;
  while(true){
    Rect<double> search_rect(q_lo - iterTime * dUnit, q_la - iterTime * dUnit, q_lo + iterTime * dUnit, q_la + iterTime * dUnit);
    num = myTree.Search(search_rect.min, search_rect.max, mySearchCallBack);
    if(totalNum + num > k){
      break;
    }
    totalNum += num;
    iterTime++;
  }
  sort(tmpTargetNode.begin(), tmpTargetNode.end(), myComp);
  vector<Node> result(tmpTargetNode.begin(), tmpTargetNode.begin() + k);
  tmpTargetNode.clear();
  return result;
}

/// 算法的最后一步:refine
/// 输入参数：
/// \param c 要查找的路径总数
/// \param rawPath IKNN中查找的所有路径
/// \param myList KNN的结果
/// \param queryPoints 存储的所有查询的点的信息
vector<pair<int, double>> refine(int c, vector<vector<Node>> rawPath, vector<vector<Node>> myList, vector<vector<double>> queryPoints){
  /// compute UB for each path
  vector<Node> tmpPath;
  vector<pair<int, double>> eachUb;
  vector<map<int, double>> simPath; //哪些点被路径遍历过并且计算过值
  //对于每条路径 计算ub 存储在eachUb中 key为路径的下标 value为ub
  for (int i = 0; i < rawPath.size(); i++){
    if(rawPath[i].empty()){
      simPath.push_back(map<int, double>());
      eachUb.push_back(pair<int, double>(0, -1));
      continue;
    }
    Node tmpNode;
    pair<int, double> tmpUb;
    map<int, double> nodes;
    double ub = 0;
    tmpPath = rawPath[i]; //对于特定的路径
    //首先 计算第一部分 match的点
    while(!tmpPath.empty()){
      tmpNode = tmpPath.back();
      tmpPath.pop_back();
      if(nodes.find(tmpNode.query_index) == nodes.end()){
        nodes[tmpNode.query_index] = calculateDist(tmpNode);
      }
      else{
        if(nodes[tmpNode.query_index] < calculateDist(tmpNode)){
          nodes[tmpNode.query_index] = calculateDist(tmpNode);
        }
      }
    }
    for(auto& x : nodes){
      ub += x.second;
    }
    simPath.push_back(nodes);
    // 其次 计算第二部分 没有match的点
    for (int j = 0; j < myList.size(); j++){
      if(nodes.find(myList[j].front().query_index) == nodes.end()){
        ub += calculateDist(myList[j].back());
      }
    }
    tmpUb.first = i;
    tmpUb.second = ub;
    eachUb.push_back(tmpUb);
  }
  /// 排序 降序
  sort(eachUb.begin(), eachUb.end(), myfunction);

  /// 开始循环寻找k个最佳路径
  vector<pair<int, double>> K_BCT;
  for (int i = 0; i < rawPath.size(); i++){//对于每一条路径
    if(rawPath[i].empty())
      continue;
    double sim = calculateSim(queryPoints, rawPath[i], simPath[i]);
    if(K_BCT.size() <= c){
      K_BCT.push_back(pair<int, double>(i, sim));
    }
    else{
      sort(K_BCT.begin(), K_BCT.end(), myfunction);
      if(sim > K_BCT.back().second){
        K_BCT.pop_back();
        K_BCT.push_back(pair<int, double>(i, sim));
      }
      if(i == rawPath.size()-1 || sim > eachUb[i + 1].second){
        return K_BCT;
      }
    }
  }
  return K_BCT;
}

/// IKNN算法主体
/// 接收参数：
/// \param c 需要选择的目标道路的条数
/// \param myTree 构造的RTREE
/// \param queryPoints
/// 返回的结果：
/// 
vector<pair<int, vector<Node>>> iK_NN(int c, RTree<Node, double, 2, double> myTree, vector<vector<double>> queryPoints){
  int k = 5;
  int deta = 2;
  //find kNN according to k and deta
  vector<vector<Node>> myList;
  while(true){
    //step 1 kNN 
    for (int i = 0; i < queryPoints.size(); i++){
      myList.push_back(qKNN(i, k, myTree, queryPoints[i][0], queryPoints[i][1]));
    }

    //step 2 merge route
    vector<vector<Node>> p(256, vector<Node>());
    vector<int> pathNum;
    int pathCount = 0;
    Node tmpNode;
    for (int i = 0; i < queryPoints.size(); i++){
      vector<Node> tmp;
      tmp = myList[i]; //tmp kNN
      while(!tmp.empty()){
        tmpNode = tmp.back();
        tmp.pop_back();
        for(int j = 0; j < tmpNode.route_index.size(); j++){
          if(p[tmpNode.route_index[j]].empty()){
            pathCount++;
            pathNum.push_back(tmpNode.route_index[j]);
          }
          p[tmpNode.route_index[j]].push_back(tmpNode); //add the node to path
        }
      }
    }

    //step 3 enough route?
    if(pathCount > c){
      //enough route 
      //then calculate LB
      vector<double> lbs(250, 0);
      priority_queue<double, vector<double>, less<double>> lb_eachPath; //最大堆
      priority_queue<double, vector<double>, greater<double>> c_Path;
      double lb = 0;
      for (int i = 0; i < pathCount; i++){
        //each queryNode
        vector<Node> tmpPath;
        tmpPath = p[pathNum[i]]; //each path
        map<int, double> m; //each query node
        for (int j = 0; j < tmpPath.size(); j++){
          if(m.find(tmpPath[j].query_index) == m.end()) //这个queryNode还没有点
            m[tmpPath[j].query_index] = calculateDist(tmpPath[j]); //only considering those points that had been matched
          else{
            if (m[tmpPath[j].query_index] < calculateDist(tmpPath[j])) //store the largest distance
              m[tmpPath[j].query_index] = calculateDist(tmpPath[j]);
          }
        }
        for(auto& x : m)
          lb += x.second;
        lb_eachPath.push(lb);
        m.clear();
      }
      for (int i = 0; i < c; i++){
        c_Path.push(lb_eachPath.top());
        lb_eachPath.pop();
      }

      //and calculate ub
      double ub = 0;
      for (int i = 0; i < myList.size(); i++){
        //printf("%lf ", calculateDist(myList[i].back()));
        ub += calculateDist(myList[i].back());
      }

      if(c_Path.top() >= ub){
        //满足条件了 需要进行refinement
        vector<pair<int, double>> K_BCT;
        K_BCT = refine(c, p, myList, queryPoints);
        vector<pair<int, vector<Node>>> returnPath;
        for (int i = 0; i < K_BCT.size(); i++){
          returnPath.push_back(pair<int, vector<Node>>(K_BCT[i].first, p[K_BCT[i].first]));
        }
        return returnPath;
      }
      else{
        k += deta;
        myList.clear();
        continue;
      }
    }
    k += deta;
    myList.clear();
  }
}


/// 主函数：
/// 1、读取数据点，构造RTREE
/// 2、随机生成50条路径，每条路径包括120个点
/// 3、随机生成查询点queryPoints
/// 4、IKNN
int main(){
  char filename[20] = "la_points.txt";
  vector<Node> points;
  vector<vector<double>> queryPoints;
  int num = 0;
  points = LoadData(filename, num);
  vector<int> randomNum;
  for (int i = 0; i < num; i++){
    randomNum.push_back(i + 1);
  }

  for (int i = 0; i < 50; i++){
    random_shuffle(randomNum.begin(), randomNum.end());
    for (int j = 0; j < 120; j ++){
      points[randomNum[j]].route_index.push_back(i);
    }
  }
  printf("Successfully Create Random Routes!\n");

  RTree<Node, double, 2, double> myTree;
  BuildRtree(points, myTree);
  printf("Successfully Build The RTree!\n");

  srand((unsigned)time(NULL));
  for (int i = 0; i < 10; i++){
    vector<double> tmp;
    tmp.push_back(rand() / double(RAND_MAX) + 33);
    tmp.push_back(rand() / double(RAND_MAX) - 119);
    queryPoints.push_back(tmp);
  }
  printf("Successfully Create Random Query Points!\n");

  /*
  RTree<Node, double, 2, double>::Iterator it;
  int itIndex = 0;
  myTree.GetFirst(it);
  while( !it.IsNull() )
  {
    Node node = *it;
    ++it;
    cout << "it[" << itIndex++ << "] " << node.lo << node.la << "\n";
  }
  */
  int c;
  printf("How many routes you wanna get:\n");
  scanf("%d", &c);
  vector<pair<int, vector<Node>>> targetPath;
  targetPath = iK_NN(c, myTree, queryPoints);
  printf("Now, we successfully finished IKNN with RTREE index.\n");
  printf("Those tajectories index are:");
  for (int i = 0; i < targetPath.size(); i++){
    printf("%d ", targetPath[i].first);
  }
  printf("\n");
  return 0;
}