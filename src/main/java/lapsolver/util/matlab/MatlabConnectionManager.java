package lapsolver.util.matlab;

import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabProxy;
import matlabcontrol.MatlabProxyFactory;

/**
 * Created by alexreinking on 6/24/14.
 */
public class MatlabConnectionManager {
    private static MatlabProxy proxy = null;

    static {
        try {
            proxy = new MatlabProxyFactory().getProxy();
        } catch (MatlabConnectionException e) {
            System.err.println("[error] Could not connect to MATLAB.");
            e.printStackTrace();
        }
    }

    private static final Object DESTRUCTOR = new Object() {
        @Override
        protected void finalize() throws Throwable {
            super.finalize();
            if (proxy != null)
                proxy.exit();
        }
    };

    public static MatlabProxy getProxy() {
        return proxy;
    }
}
