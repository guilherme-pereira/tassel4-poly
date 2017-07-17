package net.maizegenetics.gwas.NAM;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.LinkedList;
import java.util.Random;
import java.util.regex.Pattern;

public class SnpDataRandomOrder extends SnpData {
	
	static final Pattern tab = Pattern.compile("\t");
	static final byte tabByte = (byte) '\t';
	static final byte eol = (byte) '\n';
	static final byte cr = (byte) '\r';
	static LinkedList<Integer> pointerList = null;
	static FileChannel fc;
	static MappedByteBuffer mbb;
	static Random rng = new Random();
	int[] pointerArray = null;
	int currentIndex;
	
	public SnpDataRandomOrder(FileNames files) {
		super(files);
		initIndex(files.snps);
		pointerArray = getPointerArray();
		numberOfSnps = pointerArray.length;
		
		reset();
	}

	public SnpDataRandomOrder(int chromosome) {
		super(chromosome);
	}
	
	private static synchronized void initIndex(File snpFile) {
		pointerList = new LinkedList<Integer>();
		if (fc != null) {
			try {
				fc.close();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		}
		try {
			FileInputStream fis = new FileInputStream(snpFile);
			fc = fis.getChannel();
			long filesize = fc.size();
			byte[] somebytes = new byte[2048];
			mbb = fc.map(FileChannel.MapMode.READ_ONLY, 0, filesize);
			
			while (mbb.remaining() > 2048) {
				int thispos = mbb.position();
				mbb.get(somebytes);
				for (int i = 0; i < 2048; i++) {
					if (somebytes[i] == eol) {
						pointerList.add(new Integer(thispos + i + 1));
					}
				}
			}
			while (mbb.remaining() > 0) {
				if (mbb.get() == eol) pointerList.add(new Integer(mbb.position()));
			}
			
			pointerList.removeLast();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	@Override
	public boolean next() {
		currentIndex++;
		
		if (currentIndex < numberOfSnps) {
			int pointer = pointerArray[currentIndex];
			byte[] somebytes = getSomeBytes(pointer);
			int firsteol = 0;
			while (somebytes[firsteol] != cr && somebytes[firsteol] != eol) firsteol++;
			String input = new String(somebytes, 0, firsteol);
			parsedLine = tab.split(input);
			if (parsedLine.length < 15) {
				System.out.println("bad line: " + new String(somebytes));
				System.out.println("parsedLine too short at " + pointer);
			}
			return true;
		} else {
			return false;
		}
	}

	private static synchronized byte[] getSomeBytes(int pointer) {
		byte[] somebytes;
		mbb.position(pointer);
		int nleft = mbb.remaining();
		if (nleft > 500) {
			somebytes = new byte[500];
		}
		else somebytes = new byte[nleft];
		
		mbb.get(somebytes);
		return somebytes;
	}
	
	@Override
	public void reset() {
		if (pointerArray != null) {
			shuffle(pointerArray);
			currentIndex = -1; //set to -1 so that the first next sets the index to 0
		}
	}

	public static synchronized int[] getPointerArray() {
		int[] pointers = new int[pointerList.size()];
		int count = 0;
		for (Integer pointer : pointerList) pointers[count++] = pointer.intValue();
		return pointers;
	}
	
	@Override
	protected void finalize() throws Throwable {
		super.finalize();
		fc.close();
	}

	@Override
	public void findTotalSnpNumber() {
		// do nothing
	}

	@Override
	public SnpData getCopy() {
		SnpDataRandomOrder newdata = new SnpDataRandomOrder(chromosome);
		newdata.files = files;
		newdata.pointerArray = getPointerArray();
		newdata.numberOfSnps = newdata.pointerArray.length;
		return newdata;
	}

	private void shuffle(int[] array) {
	    // i is the number of items remaining to be shuffled.
		int n = array.length;
	    for (int i = n; i > 1; i--) {
	        // Pick a random element to swap with the i-th element.
	        int j = rng.nextInt(i);  // 0 <= j <= i-1 (0-based array)
	        // Swap array elements.
	        int tmp = array[j];
	        array[j] = array[i-1];
	        array[i-1] = tmp;
	    }
	}
	
}
