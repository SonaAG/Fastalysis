## Install CPU-only PyTorch
#pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# Install chatbot dependencies
#pip install streamlit transformers
#https://www.geeksforgeeks.org/gate/gate-cse-short-notes/
# run node server.js  in the backend folder first to start the API server

import streamlit as st
import requests

API_URL = "http://localhost:8000/chat"  # Change to your backend endpoint

st.title("🧬 Custom AI Chatbot")

if "messages" not in st.session_state:
    st.session_state.messages = []

user_input = st.chat_input("Ask me something about biomedical science...")

if user_input:
    st.session_state.messages.append({"role": "user", "content": user_input})

    # Send user input to your backend API
    payload = {
        "message": user_input,
        "model": st.session_state.get("selected_model", "kimi-k2")
    }
    response = requests.post(API_URL, json=payload)
    if response.ok:
        reply = response.json().get("response", "")
    else:
        reply = "Sorry, the model API is not responding."

    st.session_state.messages.append({"role": "assistant", "content": reply})

# Model selection dropdown
st.sidebar.header("Model Selection")
selected_model = st.sidebar.selectbox(
    "Choose a model",
    ["kimi-k2", "meta-llama/llama-4-scout-17b-16e-instruct", "llama-3.1-8b-instant", "deepseek-r1-distill-llama-70b", "llama-3.3-70b-versatile"],
    index=0
)
st.session_state["selected_model"] = selected_model

CHAT_STYLE = """
<style>
.chat-container {
    margin-top: 24px;
    margin-bottom: 24px;
}
.msg-row {
    display: flex;
    align-items: flex-end;
    margin-bottom: 14px;
}
.msg-row.user {
    justify-content: flex-end;
}
.msg-row.bot {
    justify-content: flex-start;
}
.avatar {
    width: 32px;
    height: 32px;
    border-radius: 50%;
    background: #1976d2;
    color: #fff;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 18px;
    font-weight: 600;
    margin: 0 10px;
    box-shadow: 0 2px 8px rgba(25, 118, 210, 0.08);
}
.bot .avatar {
    background: #7b1fa2;
}
.bubble {
    font-size: 15px;
    line-height: 1.6;
    word-break: break-word;
    border-radius: 16px;
    box-shadow: 0 2px 8px rgba(60,60,60,0.07);
    padding: 12px 18px;
}
.user .bubble {
    max-width: 65%;
    background: #e3f2fd;
    color: #222;
    border: 1px solid #e0e0e0;
    border-bottom-right-radius: 6px;
}
.bot .bubble {
    max-width: 100%;
    background: #fff;
    color: #222;
    border: none;
    border-bottom-left-radius: 6px;
    box-shadow: none;
}
</style>
"""

st.markdown(CHAT_STYLE, unsafe_allow_html=True)
st.markdown('<div class="chat-container">', unsafe_allow_html=True)

for msg in st.session_state.messages:
    if msg["role"] == "user":
        st.markdown(
            f'''
            <div class="msg-row user">
                <div class="bubble">{msg["content"]}</div>
            </div>
            ''',
            unsafe_allow_html=True
        )
    else:
        # Render bot message as Markdown inside bubble for bold/italic support
        st.markdown(
            f'''
            <div class="msg-row bot">
                <div class="bubble">''',
            unsafe_allow_html=True
        )
        st.markdown(msg["content"])
        st.markdown(
            '''</div>
            </div>
            ''',
            unsafe_allow_html=True
        )

st.markdown('</div>', unsafe_allow_html=True)
