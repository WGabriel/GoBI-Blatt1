import java.util.HashSet;

public class Test {
    int start;
    int end;

    public Test(int start, int end) {
        this.start = start;
        this.end = end;
    }


    public static void main(String[] args) {
        System.out.println("Main Started.");
        HashSet<Test> testSet1 = new HashSet<>();
        testSet1.add(new Test(2, 54));
        testSet1.add(new Test(4, 2));

        HashSet<Test> testSet2 = new HashSet<>();
        testSet2.add(new Test(4, 2));

        testSet1.retainAll(testSet2);
        for (Test t : testSet1) {
            System.out.println(t.start + " and " + t.end);
        }
    }
}
