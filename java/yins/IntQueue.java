/**
 * @file IntQueue.java
 * @author Dan Spielman <spielman@math.mit.edu>
 * @date Sun Aug  1 14:49:34 2004
 *
 * @brief integer queues of fixed maximum length
 *
 */

package yins;

/*!
 * A class implementing queues of integers
 * of fixed maximum length.
 * To be used when we will eventually use all the
 * space, so that we might as well allocate it now.
 *
 *  Relatively fast, and can be implemented quickly in C
 *
 * Warning: this class does no range checking!!!
 *
 * It is first-in first-out
*/
public class IntQueue {

    public int[] q;
    public int qLen;
    public int qPtr;

    /**
     * initialize with the maximum length
     * Note: Class does no range checking!
     *
     * @param maxLen maximum length of queue
     * @return IntQueue
     */
    public IntQueue(int maxLen) {
        q = new int[maxLen];
    }

    /**
     * Initialized the queue, and puts on entry
     *
     * @param entry first item to put on queue
     */
    public void init(int entry) {
        q[0] = entry;
        qLen = 1;
        qPtr = 0;
    }

    /**
     * Initialize the queue empty
     */
    public void init() {
        qLen = 0;
        qPtr = 0;
    }

    public boolean hasMore() {
        return (qPtr < qLen);
    }

    /**
     * Pulls off a queue item in first-in first-out order
     *
     * @return the queue item
     */
    public int pull() {
        return q[qPtr++];
    }

    /**
     * Puts an item on the queue
     *
     * @param entry entry to put on queue
     */
    public void add(int entry) {
        q[qLen++] = entry;
    }
}
