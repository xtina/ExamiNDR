package ExamiNDR.v1;

import java.io.FileNotFoundException;

public class MEME {

	public static void execute(String filepath, boolean contrastBoost,
			boolean m1, String oroot) throws FileNotFoundException {
		/*
		 * File executable = new File(mindtctPath); if(!executable.exists())
		 * throw new FileNotFoundException("Mindtct executable not found. ");
		 * 
		 * File outputDir = new File(oroot.substring(0,
		 * oroot.lastIndexOf(File.separator))); if(DEBUG)
		 * System.out.println("AbsolutePath="+outputDir.getAbsolutePath());
		 * 
		 * String[] command = null;//!contrastBoost ? {mindtctPath.trim(),
		 * arg2.trim(),
		 * {mindtctPath.trim(),arg1.trim(),arg2.trim(),filepath.trim(),arg3};
		 * List<String> argList = new Vector<String>();
		 * 
		 * command = argList.toArray(new String[argList.size()]); if(DEBUG)
		 * System.out.println("COMMAND> " +
		 * java.util.Arrays.toString(command));//"+command[0]+" "+command[1]+"
		 * "+command[2]+" "+command[3]+" "+command[4]); String[] envp = null; //
		 * should inherit the environment Runtime rt = Runtime.getRuntime();
		 */

		try {
			// Process proc = rt.exec(command, envp);
			Process proc = new ProcessBuilder("mindtct", "-m1", "-b", filepath,
					oroot).start();
			proc.waitFor();

		} catch (Exception e) {

			System.out.println("Error");
			e.printStackTrace();
		}
	}
}
