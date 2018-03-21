package plotting_v2;

import java.util.HashSet;

public class Event {
    HashSet<String> geneIds;
    int numberOfEvents;

    public Event(HashSet<String> geneIds, int numberOfEvents) {
        this.geneIds = geneIds;
        this.numberOfEvents = numberOfEvents;
    }
}
