package lapsolver.util.matlab;

import lapsolver.util.matlab.MatlabConnectionManager;
import matlabcontrol.MatlabProxy;
import org.junit.Test;
import org.junit.Assert;

public class MatlabConnectionManagerTest {
    @Test
    public void testTestConnection() throws Exception {
        MatlabProxy proxy = MatlabConnectionManager.getProxy();
        proxy.setVariable("a", 10);
        proxy.eval("a = a + 6");
        double a = ((double[]) proxy.getVariable("a"))[0];

        Assert.assertEquals(a, 16.0, Double.MIN_VALUE);
    }
}
