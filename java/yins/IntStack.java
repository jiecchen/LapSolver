/**
 * @file IntStack.java
 * @author Dan Spielman <spielman@math.mit.edu>
 * @date Tue Feb 1 2005
 *
 * @brief integer stack of fixed maximum length
 *
 */

package yins;

/*!
 * A class implementing stacks of integers
 * of fixed maximum length.
 * To be used when we will eventually use all the
 * space, so that we might as well allocate it now.
 *
 * Relatively fast, and can be implemented quickly in C
 *
 * Warning: this class does no range checking!!!
 *
 * It is last-in first-out
*/
public class IntStack {

    public int[] q;
    public int qLen;

    /**
     * initialize with the maximum length
     * Note: Class does no range checking!
     *
     * @param maxLen maximum size of stack
     * @return IntStack
     */
    public IntStack(int maxLen) {
        q = new int[maxLen];
    }

    /**
     * Initialized the stack, and puts on entry
     *
     * @param entry first item to put on stack
     */
    public void init(int entry) {
        q[0] = entry;
        qLen = 1;
    }

    /**
     * Initialize the stack empty
     *
     * @param initialize empty
     */
    public void init() {
        qLen = 0;
    }

    public boolean hasMore() {
        return (qLen > 0);
    }

    /**
     * Pulls off a stack item in last-in first-out order
     *
     * @return the stack item
     */
    public int pull() {
        return q[--qLen];
    }

    /**
     * Puts an item on the stack
     *
     * @param entry entry to put on stack
     */
    public void add(int entry) {
        q[qLen++] = entry;
    }
}
