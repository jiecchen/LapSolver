/**
 * @file   IntStack.java
 * @author Dan Spielman <spielman@math.mit.edu>
 * @date   Tue Feb 1 2005
 * 
 * @brief  integer stack of fixed maximum length
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
 * It is last-on first-off
*/
public class IntStack
{

    public int[] q;
    public int qLen;  

    /** 
     * initialize with the maximum length
     * Note: Class does no range checking!
     * 
     * @param maxLen maximum length of queue
     * 
     * @return IntStack
     */
    public IntStack(int maxLen) 
    {
	q = new int[maxLen];
    }


    /** 
     * Initialized the queue, and puts on entry
     * 
     * @param entry first item to put on queue
     */    
    public void init(int entry) 
    {
	q[0] = entry;

	qLen = 1;
	
    }

    /** 
     * Initialize the queue empty
     * 
     * @param initialize empty
     */    
    public void init() 
    {
	qLen = 0;
    }

    public boolean hasMore()
    {
	return (qLen > 0);
    }
    
    /** 
     * Pulls off a queue item in first-on first-off order
     * 
     * @return the queue item
     */
    public int pull() 
    {
	return q[--qLen];
    }

    /** 
     * Puts an item on the queue
     * 
     * @param entry entry to put on queue
     * 
     */
    public void add(int entry)
    {
	q[qLen++] = entry;
    }
	
}
