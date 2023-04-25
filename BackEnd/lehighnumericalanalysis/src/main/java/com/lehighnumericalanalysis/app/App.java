package lehighnumericalanalysis.src.main.java.com.lehighnumericalanalysis.app;

/**
 * mvn package
 * java -cp target/lehighnumericalanalysis-1.0-SNAPSHOT.jar src/main/java/com/lehighnumericalanalysis/app/App.java  *
 */
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RestController;

@SpringBootApplication
@RestController
public class App {

  @GetMapping("/api/function")
  public String function() {
    // Call the function in the .cc file here and return the result
    String result = callFunction();
    return result;
  }

  public static void main(String[] args) {
    SpringApplication.run(BackendServer.class, args);
  }
}


