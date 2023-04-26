package data;

import java.util.Collection;
import java.util.Collections;

public class Minimum_capacity implements CapacitySelector {

    @Override
    public Integer selection(Collection<Integer> sector_TV) {
        return Collections.min(sector_TV);
    }
    
}
