"""Export this task's recorded public messages, with credential redaction.

The output is a local conversation artifact. Never reads message text from
reasoning, compaction, tool, system, developer, or inter-agent records.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

THREAD = '01a06df3-c625-7a20-b6f8-a2bcb633b666'
THROUGH = "Lets zoom out, let's make a full transcript right now, lets find out what we missed, where to search, and how we find a counterexample"
TOKEN_PATTERNS = [
    re.compile(r'github(?:\\?_)?pat(?:\\?_)?[A-Za-z0-9_\\]+', re.I),
    re.compile(r'gh[pousr]_[A-Za-z0-9]{12,}'),
]


def text_content(value):
    if isinstance(value, str):
        return value
    return '\n'.join(part.get('text', '') for part in value
                     if isinstance(part, dict) and 'text' in part)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    messages, objectives = [], []
    redactions = 0
    first_timestamp = None
    cutoff = None
    seen_objectives = set()

    def redact(text):
        nonlocal redactions
        for pattern in TOKEN_PATTERNS:
            text, n = pattern.subn('[REDACTED GITHUB CREDENTIAL]', text)
            redactions += n
        return text

    with args.source.open() as stream:
        for line_number, line in enumerate(stream, 1):
            record = json.loads(line)
            payload = record.get('payload', {})
            first_timestamp = first_timestamp or record.get('timestamp')
            if record.get('type') == 'session_meta':
                assert payload.get('id') == THREAD, 'Wrong task log'
            if record.get('type') != 'event_msg':
                continue
            if payload.get('type') == 'thread_goal_updated':
                goal = payload.get('goal', {})
                objective = goal.get('objective')
                if isinstance(objective, str) and objective and objective not in seen_objectives:
                    seen_objectives.add(objective)
                    objectives.append({'timestamp': record['timestamp'],
                                       'text': redact(objective)})
                continue
            if payload.get('type') != 'item_completed':
                continue
            assert payload.get('thread_id') == THREAD, 'Cross-task message'
            item = payload.get('item', {})
            kind = item.get('type')
            if kind not in ('UserMessage', 'AgentMessage'):
                continue
            phase = item.get('phase')
            if kind == 'AgentMessage' and phase not in ('commentary', 'final_answer', 'final'):
                raise ValueError('Unrecognized public assistant phase: '+str(phase))
            original = text_content(item.get('content', []))
            messages.append({'message': len(messages)+1,
                             'timestamp_utc': record['timestamp'],
                             'role': 'user' if kind == 'UserMessage' else 'assistant',
                             'phase': phase, 'source_line': line_number,
                             'text': redact(original)})
            if kind == 'UserMessage' and original.strip() == THROUGH:
                cutoff = record['timestamp']
                break

    assert cutoff, 'Requested snapshot endpoint was not found'
    counts = Counter(message['role'] for message in messages)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    json_path = out/'CONVERSATION_TRANSCRIPT.jsonl'
    json_path.write_text(''.join(json.dumps(m, ensure_ascii=False)+'\n' for m in messages))
    md = [
        '# FASTRA H16P — full recorded visible conversation\n',
        f'Snapshot through the zoom-out request at **{cutoff}**. Times are UTC; '
        'the local date in America/Edmonton is September 4, 2026.\n',
        f'This exports every recorded public user/assistant message event in '
        f'this local task log through that request: **{len(messages)} messages** '
        f'({counts["user"]} user, {counts["assistant"]} assistant).\n',
        'Message text is preserved, except the supplied GitHub credential is '
        'redacted. Timestamps are recorded completion timestamps. Historical '
        'claims are retained and can be superseded; use the accompanying '
        'research audit for current conclusions.\n',
        'Private reasoning, system/developer instructions, automatic context, '
        'tool bodies and inter-agent messages are not part of this conversation '
        'export. Early strikes were started through the task-goal interface; '
        'the recorded task objectives are included separately below. No '
        'missing conversation with Fable, Claude or another task is invented.\n',
    ]
    for message in messages:
        phase = '' if message['role'] == 'user' else ' / '+message['phase']
        md.append(f'## M{message["message"]:03d} — {message["role"].title()}{phase} — '
                  f'{message["timestamp_utc"]}\n\n{message["text"]}\n')
    md.append('# Appendix — recorded task objectives\n\nThese are goal-interface '
              'records, not invented chat turns. Repeated copies of the same '
              'objective are listed once.\n')
    for n, objective in enumerate(objectives, 1):
        md.append(f'## Objective {n} — first recorded {objective["timestamp"]}\n\n'
                  f'{objective["text"]}\n')
    markdown_path = out/'CONVERSATION_TRANSCRIPT.md'
    markdown_path.write_text('\n'.join(md))
    for file in (json_path, markdown_path):
        assert not any(pattern.search(file.read_text()) for pattern in TOKEN_PATTERNS)
    manifest = {
        'status': 'COMPLETE_RECORDED_PUBLIC_MESSAGE_SNAPSHOT',
        'scope': 'This local task log only; no other chat history is inferred.',
        'source_log_name': args.source.name,
        'source_log_first_timestamp': first_timestamp,
        'snapshot_through_utc': cutoff,
        'exported_at_utc': datetime.now(timezone.utc).isoformat(),
        'public_messages': len(messages), 'roles': dict(counts),
        'distinct_recorded_objectives': len(objectives),
        'credential_redactions': redactions,
        'files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in (json_path, markdown_path)},
    }
    (out/'TRANSCRIPT_MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
