import streamlit as st
from datetime import datetime

def contact_page():

    st.set_page_config(page_title="Contact & Support", page_icon="📫", layout="wide")

    st.header("📫 Contact & Support")

    st.markdown("""
    If you encounter a bug, have a question, or want to request a new feature, please open an issue on our GitHub repository.  
    
    This helps us track and resolve problems transparently.

    You can also include screenshots or error messages in your issue.
    """)

    # --- Replace this with your actual repo URL ---
    GITHUB_REPO = "https://github.com/YourUsername/YourRepoName"

    st.link_button("🛠️  Open an issue on GitHub", f"{GITHUB_REPO}/issues/new/choose")

    st.divider()

    st.subheader("💬 Other ways to reach out")
    st.markdown("""
    - **Email:** [support@example.org](mailto:support@example.org)  
    - **Documentation:** [Project Wiki](https://github.com/YourUsername/YourRepoName/wiki)
    """)

    st.caption(f"© {datetime.now().year} Daniel Paysan — built with ❤️ and Streamlit")
