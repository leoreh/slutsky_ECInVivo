import sys

pdf_path = sys.argv[1]
txt_path = sys.argv[2]
try:
    import fitz
    text = ""
    with fitz.open(pdf_path) as doc:
        for page in doc:
            text += page.get_text()
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write(text)
    print("PyMuPDF success")
except ImportError:
    try:
        import PyPDF2
        with open(pdf_path, 'rb') as f:
            reader = PyPDF2.PdfReader(f)
            text = '\n'.join([p.extract_text() for p in reader.pages])
        with open(txt_path, "w", encoding="utf-8") as out:
            out.write(text)
        print("PyPDF2 success")
    except Exception as e:
        print(f"Error: {e}")
