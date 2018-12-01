#pragma once

#include <thread>
using namespace std::placeholders;

MTS_NAMESPACE_BEGIN

/**
 * The block scheduler is in charge of dispatching a block to each thread for processing. 
 * 
 * It takes a user function that defines how to process each block and attempts to execute 
 * it in parallel until no more blocks are available. 
 */


/************************************************************************
 * Block Scheduler
 ************************************************************************/
class BlockScheduler {

  /************************************************************************
   * Local Class - Block Thread
   ************************************************************************/
public:

  class BlockThread : public Thread {
  public:
    //Store the lambda
    typedef std::function<void(int, int)> ComputeBlockFunction;
    ComputeBlockFunction localF;
    int localTID = -1;
    BlockScheduler &localParent;

    //Constructor
    BlockThread(const std::string &threadName, ComputeBlockFunction f, int tid, BlockScheduler &parent)
        : Thread(threadName), localParent(parent) {
      this->localF = f;
      this->localTID = tid;
      this->setPriority(EThreadPriority::ENormalPriority);
    }

    virtual void run() {
      while (true) {
        std::tuple<int, int> blockIdx = localParent.getBlockIdx();
        int begIdx = std::get<0>(blockIdx);
        int endIdx = std::get<1>(blockIdx);
        if (begIdx < 0)
          break;
        for (int id = begIdx; id < endIdx; id++) {
          localF(id, localTID);
        }
      }
    }
  };


  /************************************************************************
   * Constructor
   ************************************************************************/
public:

  BlockScheduler(int numBlocks, int numThreads, int granuarity = 1) :
      numBlocks(numBlocks),
      numThreads(numThreads),
      blockIdx(0),
      granuarity(granuarity) {
    mutex = new Mutex();
  }


  /************************************************************************
   * Typedefs - ComputeBlockFunction
   ************************************************************************/
public:

  typedef std::function<void(int, int)> ComputeBlockFunction;


  /************************************************************************
   * Public Functions
   ************************************************************************/
public:

  /**
   * Runs a ComputeBlockFunctino for numBlocks on numThreads
   */
  void run(ComputeBlockFunction f) {
    ref_vector<BlockThread> group;
    for (int tid = 0; tid < numThreads; ++tid) {
      ref<BlockThread> bThread = new BlockThread("BLOCKTHREAD" + std::to_string(tid), f, tid, *this);
      group.push_back(bThread);
      bThread->start();
    }

    for (int tid = 0; tid < numThreads; ++tid) {
      group[tid]->join();
    }
  }

  /**
   * Return a unique block for each thread.
   * Return a negative number when no blocks are available.
   */
  std::tuple<int, int> getBlockIdx() {
    LockGuard lock(mutex);
    if (blockIdx >= numBlocks)
      return std::make_tuple(-1, -1);

    int v = blockIdx;
    int v2 = std::min(blockIdx + granuarity, numBlocks);
    blockIdx = v2;
    return std::make_tuple(v, v2);
  }


  /************************************************************************
   * Destructor
   ************************************************************************/
public:

  ~BlockScheduler() {

  }


  /************************************************************************
   * Private Class Variables
   ************************************************************************/
private:

  int numBlocks;
  int numThreads;
  int blockIdx;
  int granuarity;
  ref <Mutex> mutex;
};

MTS_NAMESPACE_END
